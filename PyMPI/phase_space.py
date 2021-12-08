from mpi4py import MPI
import math
import pandas as pd
import numpy as np
import integration
import constants as k


def evolution(s, dist):
    postTime = dist / 164.35
    s["b"] += postTime
    s["bT"] += postTime

    if s["chirp"] > 0:
        s["VzDist"] = math.sqrt(
            1 / ((1 / pow(s["hHeight"], 2)) + pow((s["b"] / s["zDist"]), 2)))

    if s["chirpT"] > 0:
        s["VxDist"] = math.sqrt(
            1 / ((1 / pow(s["hDepthVel"], 2)) + pow((s["bT"] / s["xDist"]), 2)))

    if s["chirp"] < 0:
        s["zDist"] = s["b"] / \
            (math.sqrt((1 / pow(s["VzDist"], 2)) - (1 / pow(s["hHeight"], 2))))

    if s["chirpT"] < 0:
        s["xDist"] = s["bT"] / \
            (math.sqrt((1 / pow(s["VxDist"], 2)) -
             (1 / pow(s["hDepthVel"], 2))))

    s["chirp"] = s["b"] * pow(s["VzDist"] / s["zDist"], 2)
    s["chirpT"] = s["bT"] * pow(s["VxDist"] / s["xDist"], 2)
    s["hWidth"] = math.sqrt(
        1 / ((1 / pow(s["zDist"], 2)) - pow(s["chirp"] / s["VzDist"], 2)))
    s["hDepth"] = math.sqrt(
        1 / ((1 / pow(s["xDist"], 2)) - pow(s["chirpT"] / s["VxDist"], 2)))

    s["zC"] += s["VzC"] * postTime
    s["xC"] += s["hDepthVel"] * postTime


def split(s, startSplit, numSplits, pulseEnergy, pool):
    splitSpaces = np.array(numSplits)
    #pool.map(square, range(1000))
    args = [(s, i+1)
            for i in range(startSplit, startSplit+numSplits)]
    #print("args", args)
    intensityMultipliers = np.asarray(
        pool.starmap(integration.get_intensity, args))
    # intensityMultipliers = np.asarray(
    #    [integration.get_intensity(s, numSplits, i+1) for i in range(startSplit, startSplit + numSplits)])
    splitHHeight = 0
    splitVzDist = 0
    splitB = 0
    splitZDist = 0
    splitChirp = 0
    splitHHeight = s["hHeight"] / k.SPLIT_NUM
    splitVzDist = s["VzDist"] / k.SPLIT_NUM

    splitB = s["zDist"] * \
        math.sqrt((1 / pow(splitVzDist, 2)) - (1 / pow(splitHHeight, 2)))
    splitZDist = splitB / \
        (math.sqrt((1 / pow(splitVzDist, 2)) - (1 / pow(splitHHeight, 2))))
    splitChirp = splitVzDist * \
        math.sqrt((1 / pow(s["zDist"], 2)) - (1 / pow(s["hWidth"], 2)))

    splitSpacesList = [[s["hWidth"], splitHHeight, splitVzDist, splitZDist, splitChirp, splitB, pulseEnergy*intensityMultipliers[j-startSplit], intensityMultipliers[j-startSplit], s["hDepth"], s["hDepthVel"], s["VxDist"],
                        s["xDist"], s["chirpT"], s["bT"], s["hHeight"]-(s["hHeight"]*2/k.SPLIT_NUM)*(j+0.5), (s["hHeight"] - (s["hHeight"] * 2 / k.SPLIT_NUM)*(j + 0.5)) / s["chirp"], s["VxC"], s["xC"]] for j in range(startSplit, startSplit + numSplits)]

    splitSpaces = pd.DataFrame(splitSpacesList, columns=k.COLUMNS)
    return splitSpaces


def evolutionWithInteraction(sharedMem, pulses, startSplit, numSplits, splitSize, time):
    # reiniitalize MPI to reduce method arguments
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # assign space values for NN to sharedMemory
    for i in range(numSplits):
        sharedMem[(startSplit+i)*3] = pulses["pulseEnergy"].iloc[i]
        sharedMem[(startSplit+i)*3+1] = pulses["VzC"].iloc[i]
        sharedMem[(startSplit+i)*3+2] = pulses["zC"].iloc[i]
    comm.Barrier()

    # create window for one-way communication (other phase-space parameters)
    oneWayBuf = pulses.to_numpy()
    oneWayBufRecv = np.zeros(18*numSplits, dtype='d')
    oneWayWin = MPI.Win.Create(oneWayBuf, 1, MPI.INFO_NULL, MPI.COMM_WORLD)
    oneWayWin.Lock(rank)
    oneWayWin.Put(oneWayBuf, target_rank=rank)
    oneWayWin.Unlock(rank)

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

    comm.Barrier()
    print(neighbors)

    # Retrive neighbor's values
    neighborSpaces = {}
    for key, value in neighbors.items():
        for val in value:
            neighborSpaces[val] = []

    # Find corresponding ranks for neighbors
    neighborRanks = {}
    for key in neighborSpaces.keys():
        targetRank = key//splitSize
        if targetRank in neighborRanks:
            neighborRanks[targetRank].append(key)
        else:
            neighborRanks[targetRank] = [key]

    # Retrieve neighbors from ranks
    for key, val in neighborRanks.items():
        oneWayWin.Lock(key)
        oneWayWin.Get(oneWayBufRecv, target_rank=key)
        oneWayWin.Unlock(key)
        spaceMatrix = oneWayBufRecv.reshape((-1, 18))
        for ind in val:
            neighborSpaces[ind] = spaceMatrix[ind % splitSize, :]

    # Next step: calculate force effect on location
    forces = [np.sum([integration.getForce(
        pd.Series(j, index=k.COLUMNS), pulses.iloc[i-startSplit]["zC"], pulses.iloc[i-startSplit]["VzC"], pulses.iloc[i-startSplit]["pulseEnergy"]) for j in neighbors[i]], axis=0) for i in range(startSplit, startSplit+numSplits)]

    forces = pd.DataFrame([force if type(force) is np.ndarray else np.array([0.0, 0.0])
                           for force in forces], columns=["zC", "VzC"])
    if rank == 1:
        print(pulses[["zC", "VzC"]])

    pulses["zC"] += forces["zC"]*time/1E15
    pulses["VzC"] += forces["VzC"]*time/1E15

    if rank == 1:
        print(pulses[["zC", "VzC"]])
