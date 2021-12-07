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
