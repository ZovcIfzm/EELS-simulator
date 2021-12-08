
import multiprocessing
import numpy as np

import math
import constants as k
import decimal

# takes in base (singular) phase space


def drange(x, y, jump):
    while x < y:
        yield float(x)
        x += jump


def get_intensity(s, sectionNum):
    # Gets intensity % proportionally to 1 (like if its gets .5 its 50 % of total intensity)
    # search with xSearch & ySearch = +- 5.803*hWidth or hHeight to get the total intensity of the phase space(equal to 1)
    ySearchLB = -k.CATCH_FACTOR * s["hHeight"] + \
        ((k.CATCH_FACTOR * s["hHeight"] *
         2.0 / k.SPLIT_NUM) * (sectionNum - 1))
    ySearchUB = k.CATCH_FACTOR * s["hHeight"] - \
        (k.CATCH_FACTOR * s["hHeight"] * 2.0 /
         k.SPLIT_NUM) * (k.SPLIT_NUM - sectionNum)
    xSearchLB = -k.CATCH_FACTOR * s["hWidth"]
    xSearchUB = k.CATCH_FACTOR * s["hWidth"]

    # Convert to intensity_integration parameters
    xHalfRange = (xSearchUB - xSearchLB) / 2
    yHalfRange = (ySearchUB - ySearchLB) / 2

    xOffset = (xSearchUB + xSearchLB) / 2
    yOffset = (ySearchUB + ySearchLB) / 2

    return intensity_integration(s, xHalfRange, yHalfRange, xOffset, yOffset)


def intensity_integration(s, xHalfRange, yHalfRange, xOffset, yOffset):
    accuracyY = 2 * yHalfRange / 199
    accuracyX = 2 * xHalfRange / 199
    x = -xHalfRange + xOffset + accuracyX / 2
    y = -yHalfRange + yOffset + accuracyY / 2

    intensityValue = 0

    while (y < yHalfRange + yOffset):
        while (x < xHalfRange + xOffset):
            intensityValue += accuracyX * accuracyY * intensity(s, x, y)
            x += accuracyX
        y += accuracyY
        x = -xHalfRange + xOffset + accuracyX / 2
    return intensityValue


def intensity(s, x, y):
    negTwohWidthsq = -2 * s["hWidth"] * s["hWidth"]
    twoVzIntDistsq = 2 * s["VzDist"] * s["VzDist"]
    twoPIhWidthVzIntDist = 2 * math.pi * (s["hWidth"] * s["VzDist"])
    return math.exp((x * x / (negTwohWidthsq)) - ((y - s["chirp"] * x) * (y - s["chirp"] * x) / (twoVzIntDistsq))) / (twoPIhWidthVzIntDist)


def force(q1, q2, x1, y1, x2, y2):
    # 1 is particle, 2 is other dist
    f = k.KC*q1*q2/(math.pow(x1-x2, 2)+math.pow(y1-y2, 2))
    angle = math.atan((x2-x1)/(y2-y1))
    return [math.sin(angle)*f, math.cos(angle)*f]


def getForce(s, refX, refY, refEnergy):
    ySearchLB = -k.CATCH_FACTOR * s["hHeight"]
    ySearchUB = k.CATCH_FACTOR * s["hHeight"]
    xSearchLB = -k.CATCH_FACTOR * s["hWidth"]
    xSearchUB = k.CATCH_FACTOR * s["hWidth"]

    # Convert to intensity_integration parameters
    xHalfRange = (xSearchUB - xSearchLB) / 2
    yHalfRange = (ySearchUB - ySearchLB) / 2

    xOffset = (xSearchUB + xSearchLB) / 2
    yOffset = (ySearchUB + ySearchLB) / 2

    accuracyY = 2 * yHalfRange / 199
    accuracyX = 2 * xHalfRange / 199
    x = -xHalfRange + xOffset + accuracyX / 2
    y = -yHalfRange + yOffset + accuracyY / 2

    intensities = np.asarray([[force(refEnergy, accuracyX*accuracyY*intensity(s, i, j), refX, refY, s["zC"], s["VzC"]) for j in list(drange(
        x, xHalfRange+xOffset, accuracyX))] for i in list(drange(y, yHalfRange + yOffset, accuracyY))])

    netForce = np.sum(np.sum(intensities, axis=0), axis=0)
    return netForce


def pixelSum(s):
    pixelArray = np.zeros((k.NUM_PIXELS))
    lowestXC = s.iloc[len(s)-1]["xC"]
    highestXC = s.iloc[0]["xC"]
    xCDist = abs(highestXC - lowestXC) / len(pixelArray)

    pixelArray = np.asarray([pixelSumHelper(s, p, lowestXC, xCDist)
                            for p in range(len(pixelArray))])
    pixelArray = pixelArray / np.sqrt(np.sum(pixelArray**2))

    return pixelArray


def pixelSumHelper(s, p, lowestXC, xCDist):
    print('pixel', p)
    returnVal = 0
    for i in range(len(s)):
        if (lowestXC + p * xCDist < s.iloc[i]["xC"] + 5 * s.iloc[i]["hDepth"] and s.iloc[i]["xC"] - 5 * s.iloc[i]["hDepth"] < lowestXC + (p + 1.0) * xCDist):
            returnVal += x_integration(s.iloc[i], lowestXC + p * xCDist, lowestXC + (
                p + 1.0) * xCDist) * s.iloc[i]["intensityMultiplier"]

    return p, returnVal


def x_integration(s, xLeftLim, xRightLim):
    # was 19 not 5
    accuracyY = 2 * s["hHeight"] / 5
    accuracyX = (xRightLim - xLeftLim) / 5
    x = xLeftLim
    y = -s["hHeight"] + accuracyY / 2
    intensityValue = 0
    while (y < s["hHeight"]):
        while (x < xRightLim):
            intensityValue += accuracyX * accuracyY * \
                intensity(s, x - s["xC"], y)
            x += accuracyX
        y += accuracyY
        x = xLeftLim

    return intensityValue
