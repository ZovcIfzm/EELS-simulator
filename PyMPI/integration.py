
import multiprocessing

import math
import constants as k

# takes in base (singular) phase space


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
