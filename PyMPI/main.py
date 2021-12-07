from mpi4py import MPI
import multiprocessing

import math
import numpy as np
import pandas as pd
from hilbertcurve.hilbertcurve import HilbertCurve

import phase_space as ps
import constants as k

# initialize mpi
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# set up multiprocessing
try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 2   # arbitrary default
print("num cpus:", cpus)

pool = multiprocessing.Pool(processes=cpus)

# create shared memory for initial split
itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    nbytes = k.SPLIT_NUM*(itemsize)*3
else:
    nbytes = 0

win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)

buf, itemsize = win.Shared_query(0)
assert itemsize == MPI.DOUBLE.Get_size()
sharedMem = np.ndarray(buffer=buf, dtype='d', shape=(3*k.SPLIT_NUM,))

# split the phase space... unless more cores than splits
# TODO, create the else statement... which starts after shatters I guess?
if rank < k.SPLIT_NUM:
    splitSize = max(math.floor(k.SPLIT_NUM/size), 1)
    startSplit = rank*splitSize
    numSplits = splitSize
    if startSplit + splitSize < k.SPLIT_NUM and startSplit + 2*splitSize > k.SPLIT_NUM:
        numSplits = k.SPLIT_NUM-startSplit

    initPulse = k.INIT_PS
    pulses = ps.split(initPulse, startSplit, numSplits, k.PULSE_ENERGY, pool)
    # assign space values for NN to sharedMemory
    for i in range(numSplits):
        sharedMem[(startSplit+i)*3] = pulses["pulseEnergy"].iloc[i]
        sharedMem[(startSplit+i)*3+1] = pulses["VzC"].iloc[i]
        sharedMem[(startSplit+i)*3+2] = pulses["zC"].iloc[i]
    comm.Barrier()

    # find nearest neighbors
    # convert values to dataframe, run parallel operations
    df = pd.DataFrame(sharedMem.reshape((k.SPLIT_NUM, 3)),
                      columns=["pulseEnergy", "VzC", "zC"])

    neighbors = {}

    for i in range(startSplit, startSplit+numSplits):
        limitDf = pd.DataFrame(df["pulseEnergy"]*df["pulseEnergy"].iloc[i] /
                               (pow(df["VzC"]-df["VzC"].iloc[i], 2) +
                                pow(df["zC"]-df["zC"].iloc[i], 2)), columns=["limit"])
        limitDf = limitDf[~limitDf.isin([np.nan, np.inf, -np.inf]).any(1)]
        limitDf = limitDf[limitDf["limit"] > k.NN_LIMIT]
        neighbors[i] = limitDf.index.tolist()

    print(neighbors)

    # map with space-filling curve
    p, n = 10, 2
    hilbertCurve = HilbertCurve(p, n)
    pointsDf = df[["VzC", "zC"]]
    if len(df.index) != 1:
        normPointsDf = (pointsDf-pointsDf.min()) / \
            (pointsDf.max()-pointsDf.min())
        intPointsDf = (normPointsDf*1E3).astype(int)

        points = intPointsDf.to_numpy()

        distances = hilbertCurve.distances_from_points(points)
        sortIndices = np.argsort(distances)

        print(distances)

'''
loop:
    Max over spaces-can use edges? to find bucket dimensions
    shared memory
    (first) nearest neighbors, parallel
        (?) Create buckets for nearest neighbor potentials.such that potentials for a bucket are adjacent buckets.
        (all broadcast, or send to single machine. single machine defines buckets, and bucket potentials, sends to appropriate processors)
    Create space filling curve, based on number of neighbors as weight.
    Divide among cpus the phase spaces & their potential buckets
    Each phase space calculates which potential bucket to take into account.
    (first) calculate intensity based effect (later) Calculate integrated effect 
    Calculate natural phase space changes
shatter
repeat above loop
reduce on 1D
'''
'''

data = None

if rank == 0:
    data = {'a': 7, 'b': 3.14}
    comm.send(data, dest=1, tag=11)
elif rank == 1:
    data = comm.recv(source=0, tag=11)
'''
