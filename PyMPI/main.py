from mpi4py import MPI
import multiprocessing

import math
import numpy as np
import pandas as pd
from hilbertcurve.hilbertcurve import HilbertCurve

import integration
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

# read spectrum
spectrum = pd.read_csv("spectrum.csv")
#print("spec:", spectrum)

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

    ps.evolutionWithInteraction(
        sharedMem, pulses, startSplit, numSplits, splitSize, 500)

    '''
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
