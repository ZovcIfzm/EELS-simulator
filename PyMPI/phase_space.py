from mpi4py import MPI
import math
import pandas as pd
import numpy as np
import integration
import constants as k
import copy
from hilbertcurve.hilbertcurve import HilbertCurve


def evolution(s, dist):
    postTime = dist / 164.35
    s["b"] += postTime
    s["bT"] += postTime

    if s["chirp"].iloc[0] > 0:
        s["VzDist"] = (1 / ((1 / pow(s["hHeight"], 2)) +
                       pow((s["b"] / s["zDist"]), 2)))**0.5

    if s["chirpT"].iloc[0] > 0:
        s["VxDist"] = (1 / ((1 / pow(s["hDepthVel"], 2)) +
                       pow((s["bT"] / s["xDist"]), 2)))**0.5

    if s["chirp"].iloc[0] < 0:
        s["zDist"] = s["b"] / \
            ((1 / pow(s["VzDist"], 2)) - (1 / pow(s["hHeight"], 2)))**0.5

    if s["chirpT"].iloc[0] < 0:
        s["xDist"] = s["bT"] / \
            ((1 / pow(s["VxDist"], 2)) -
             (1 / pow(s["hDepthVel"], 2)))**0.5

    s["chirp"] = s["b"] * pow(s["VzDist"] / s["zDist"], 2)
    s["chirpT"] = s["bT"] * pow(s["VxDist"] / s["xDist"], 2)
    s["hWidth"] = (1 / ((1 / pow(s["zDist"], 2)) -
                   pow(s["chirp"] / s["VzDist"], 2)))**0.5
    s["hDepth"] = (1 / ((1 / pow(s["xDist"], 2)) -
                   pow(s["chirpT"] / s["VxDist"], 2)))**0.5

    s["zC"] += s["VzC"] * postTime
    s["xC"] += s["hDepthVel"] * postTime


def split(s, startSplit, numSplits, pulseEnergy, pool):
    splitSpaces = np.array(numSplits)
    # pool.map(square, range(1000))
    args = [(s, i+1)
            for i in range(startSplit, startSplit+numSplits)]
    # print("args", args)
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


def shatter(s, spectrum):
    spaces = len(spectrum)
    print("spaces", spaces)
    print("s length", len(s))
    newVzDist = s["VzDist"] / spaces
    newChirp = newVzDist * \
        ((1 / pow(s["zDist"], 2)) - (1 / pow(s["hWidth"], 2)))**0.5
    newB = newChirp * pow(s["zDist"] / newVzDist, 2)
    shatteredPulses = pd.concat([pd.concat([s["hWidth"], s["hHeight"] / spaces, newVzDist, s["zDist"], newChirp, newB, s["pulseEnergy"] * spectrum.iloc[i, 1] / k.BASE_TOTAL, s["intensityMultiplier"] * spectrum.iloc[i, 1],
                                            s["hDepth"], s["hDepthVel"], s["VxDist"], s["xDist"], s["chirpT"], s["bT"], s["hHeight"] + (spectrum.iloc[i, 0] / 1117), s["zC"], s["VxC"], s["xC"]], axis=1, keys=k.COLUMNS) for i in range(spaces)], ignore_index=True)

    return shatteredPulses


def evolutionWithInteraction(sharedMem, pulses, startSplit, numSplits, splitSize, dist, perfect=False):
    totTime = dist / 164.35
    passedTime = 0.0
    timeInc = min(1, totTime)
    while passedTime < totTime:
        passedTime += timeInc
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

        if size != 1:
            # create window for one-way communication (other phase-space parameters)
            oneWayBuf = pulses.to_numpy()
            oneWayBufRecv = np.zeros(18*numSplits, dtype='d')
            oneWayWin = MPI.Win.Create(
                oneWayBuf, 1, MPI.INFO_NULL, MPI.COMM_WORLD)
            oneWayWin.Lock(rank)
            oneWayWin.Put(oneWayBuf, target_rank=rank)
            oneWayWin.Unlock(rank)

            # find nearest neighbors
            # convert values to dataframe, run parallel operations
            df = pd.DataFrame(sharedMem.reshape((k.SPLIT_NUM, 3)),
                              columns=["pulseEnergy", "VzC", "zC"])
        else:
            df = pulses

        neighbors = {}

        if perfect:
            for i in range(startSplit, startSplit+numSplits):
                limitDf = pd.DataFrame(df["pulseEnergy"]*df["pulseEnergy"].iloc[i] /
                                       (pow(df["VzC"]-df["VzC"].iloc[i], 2) +
                                        pow(df["zC"]-df["zC"].iloc[i], 2)), columns=["limit"])
                limitDf = limitDf[~limitDf.isin(
                    [np.nan, np.inf, -np.inf]).any(1)]
                print("df", limitDf)
                print("dflim", k.NN_LIMIT)
                limitDf = limitDf[limitDf["limit"] > k.NN_LIMIT]
                neighbors[i] = limitDf.index.tolist()
        else:
            # calculate hilbert curve
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
                sortDist = sorted(distances)

                for i in range(startSplit, startSplit+numSplits):
                    distance = distances[i]
                    index = np.where(sortIndices == i)[0][0]
                    j = 1
                    neighbors[i] = []
                    while index - j >= 0 and abs(sortDist[index-j]-distance) < k.HILBERT_THRESHOLD:
                        neighbors[i].append(sortIndices[index-j])
                        j = j + 1
                    j = 1
                    while index + j < len(df) and abs(sortDist[index+j]-distance) < k.HILBERT_THRESHOLD:
                        neighbors[i].append(sortIndices[index+j])
                        j = j + 1

        comm.Barrier()

        if size != 1:
            # Retrive neighbor's values
            neighborSpaces = {}
            for key, value in neighbors.items():
                for val in value:
                    neighborSpaces[val] = []

            # Find corresponding ranks for neighbors
            neighborRanks = {}
            for key in neighborSpaces.keys():
                targetRank = min(key//splitSize, size-1)
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
                pd.Series(neighborSpaces[j], index=k.COLUMNS), pulses.iloc[i-startSplit]["zC"], pulses.iloc[i-startSplit]["VzC"], pulses.iloc[i-startSplit]["pulseEnergy"]) for j in neighbors[i]], axis=0) for i in range(startSplit, startSplit+numSplits)]

        else:
            forces = [np.sum([integration.getForce(
                pd.Series(df.iloc[j], index=k.COLUMNS), pulses.iloc[i-startSplit]["zC"], pulses.iloc[i-startSplit]["VzC"], pulses.iloc[i-startSplit]["pulseEnergy"]) for j in neighbors[i]], axis=0) for i in range(startSplit, startSplit+numSplits)]

        forces = pd.DataFrame([force if type(force) is np.ndarray else np.array([0.0, 0.0])
                               for force in forces], columns=["zC", "VzC"])

        pulses["zC"] += forces["zC"]*timeInc/1E15
        pulses["VzC"] += forces["VzC"]*timeInc/1E15
        evolution(pulses, timeInc*164.35)


def evolutionWithNoSharing(pulses, startSplit, numSplits, splitSize, dist, perfect=False):
    totTime = dist / 164.35
    passedTime = 0.0
    timeInc = min(1, totTime)
    while passedTime < totTime:
        passedTime += timeInc
        # reiniitalize MPI to reduce method arguments
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # create window for one-way communication (other phase-space parameters)

        if size != 1:
            oneWayBuf = pulses.to_numpy().reshape((-1, 1))
            oneWayBufRecv = np.zeros(18*numSplits, dtype='d')
            oneWayWin = MPI.Win.Create(
                oneWayBuf, 1, MPI.INFO_NULL, MPI.COMM_WORLD)
            oneWayWin.Lock(rank)
            oneWayWin.Put(oneWayBuf, target_rank=rank)
            oneWayWin.Unlock(rank)

        # find nearest neighbors
        # retrieve all neighbor's values
        # convert values to dataframe, run parallel operations

        if size != 1:
            allPulses = []
            for i in range(size):
                oneWayWin.Lock(i)
                oneWayWin.Get(oneWayBufRecv, target_rank=i)
                oneWayWin.Unlock(i)
                allPulses.append(copy.deepcopy(oneWayBufRecv))

            allPulses = np.asarray(allPulses).reshape((-1, 18))

            df = pd.DataFrame(allPulses,
                              columns=k.COLUMNS)
        else:
            df = pulses

        neighbors = {}
        if perfect:
            for i in range(startSplit, startSplit+numSplits):
                limitDf = pd.DataFrame(df["pulseEnergy"]*df["pulseEnergy"].iloc[i] /
                                       (pow(df["VzC"]-df["VzC"].iloc[i], 2) +
                                        pow(df["zC"]-df["zC"].iloc[i], 2)), columns=["limit"])
                limitDf = limitDf[~limitDf.isin(
                    [np.nan, np.inf, -np.inf]).any(1)]
                limitDf = limitDf[limitDf["limit"] > k.NN_LIMIT]
                neighbors[i] = limitDf.index.tolist()
        else:
            # calculate hilbert curve
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
                sortDist = sorted(distances)

                for i in range(startSplit, startSplit+numSplits):
                    distance = distances[i]
                    index = np.where(sortIndices == i)[0][0]
                    j = 1
                    neighbors[i] = []
                    while index - j >= 0 and abs(sortDist[index-j]-distance) < k.HILBERT_THRESHOLD:
                        neighbors[i].append(sortIndices[index-j])
                        j = j + 1
                    j = 1
                    while index + j < len(df) and abs(sortDist[index+j]-distance) < k.HILBERT_THRESHOLD:
                        neighbors[i].append(sortIndices[index+j])
                        j = j + 1

        comm.Barrier()
        # print(neighbors)

        # Retrive neighbor's values
        neighborSpaces = {}
        for key, value in neighbors.items():
            for val in value:
                neighborSpaces[val] = []

        # Find corresponding ranks for neighbors
        neighborRanks = {}
        for key in neighborSpaces.keys():
            targetRank = min(key//splitSize, size-1)
            if targetRank in neighborRanks:
                neighborRanks[targetRank].append(key)
            else:
                neighborRanks[targetRank] = [key]

        print(neighbors)
        if size != 1:
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
                pd.Series(neighborSpaces[j], index=k.COLUMNS), pulses.iloc[i-startSplit]["zC"], pulses.iloc[i-startSplit]["VzC"], pulses.iloc[i-startSplit]["pulseEnergy"]) for j in neighbors[i]], axis=0) for i in range(startSplit, startSplit+numSplits)]
        else:
            forces = [np.sum([integration.getForce(
                pd.Series(pulses.iloc[j], index=k.COLUMNS), pulses.iloc[i-startSplit]["zC"], pulses.iloc[i-startSplit]["VzC"], pulses.iloc[i-startSplit]["pulseEnergy"]) for j in neighbors[i]], axis=0) for i in range(startSplit, startSplit+numSplits)]
        forces = pd.DataFrame([force if type(force) is np.ndarray else np.array([0.0, 0.0])
                               for force in forces], columns=["zC", "VzC"])

        pulses["zC"] += forces["zC"]*timeInc/1E15
        pulses["VzC"] += forces["VzC"]*timeInc/1E15
        evolution(pulses, timeInc*164.35)


def evolutionWithFullInteraction(fullSharedMem, pulses, startSplit, numSplits, splitSize, dist, perfect=False):
    totTime = dist / 164.35
    passedTime = 0.0
    timeInc = min(1, totTime)
    while passedTime < totTime:
        passedTime += timeInc
        # reiniitalize MPI to reduce method arguments
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        # assign space values for NN to sharedMemory
        for i in range(numSplits):
            pulse = pulses.iloc[i].to_numpy()
            for j in range(len(pulse)):
                fullSharedMem[(startSplit+i)*18+j] = pulse[j]

        comm.Barrier()

        # find nearest neighbors
        # convert values to dataframe, run parallel operations
        df = pd.DataFrame(fullSharedMem.reshape((k.SPLIT_NUM, 18)),
                          columns=k.COLUMNS)

        neighbors = {}
        if perfect:
            for i in range(startSplit, startSplit+numSplits):
                limitDf = pd.DataFrame(df["pulseEnergy"]*df["pulseEnergy"].iloc[i] /
                                       (pow(df["VzC"]-df["VzC"].iloc[i], 2) +
                                        pow(df["zC"]-df["zC"].iloc[i], 2)), columns=["limit"])
                limitDf = limitDf[~limitDf.isin(
                    [np.nan, np.inf, -np.inf]).any(1)]
                limitDf = limitDf[limitDf["limit"] > k.NN_LIMIT]
                neighbors[i] = limitDf.index.tolist()
        else:
            # calculate hilbert curve
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
                sortDist = sorted(distances)

                for i in range(startSplit, startSplit+numSplits):
                    distance = distances[i]
                    index = np.where(sortIndices == i)[0][0]
                    j = 1
                    neighbors[i] = []
                    while index - j >= 0 and abs(sortDist[index-j]-distance) < k.HILBERT_THRESHOLD:
                        neighbors[i].append(sortIndices[index-j])
                        j = j + 1
                    j = 1
                    while index + j < len(df) and abs(sortDist[index+j]-distance) < k.HILBERT_THRESHOLD:
                        neighbors[i].append(sortIndices[index+j])
                        j = j + 1

        comm.Barrier()
        print(neighbors)
        # Next step: calculate force effect on location
        forces = [np.sum([integration.getForce(
            pd.Series(df.iloc[j], index=k.COLUMNS), pulses.iloc[i-startSplit]["zC"], pulses.iloc[i-startSplit]["VzC"], pulses.iloc[i-startSplit]["pulseEnergy"]) for j in neighbors[i]], axis=0) for i in range(startSplit, startSplit+numSplits)]

        forces = pd.DataFrame([force if type(force) is np.ndarray else np.array([0.0, 0.0])
                               for force in forces], columns=["zC", "VzC"])

        pulses["zC"] += forces["zC"]*timeInc/1E15
        pulses["VzC"] += forces["VzC"]*timeInc/1E15
        evolution(pulses, timeInc*164.35)


def magLens(s, power):
    s["chirpT"] -= pow(power, 2) * k.MAG_LENS_COEFFICIENT
    s["xDist"] = (1 / ((1 / pow(s["hDepth"], 2)) +
                  pow(s["chirpT"] / s["VxDist"], 2)))**0.5
    s["bT"] = s["chirpT"] * pow(s["xDist"] / s["VxDist"], 2)
    s["hDepthVel"] = (1 / ((1 / pow(s["VxDist"], 2)) -
                      pow(s["bT"] / s["xDist"], 2)))**0.5
