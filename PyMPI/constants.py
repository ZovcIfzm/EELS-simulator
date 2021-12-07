import pandas as pd
import numpy as np

SPLIT_NUM = 4
NUM_CORES = 36
CATCH_FACTOR = 4
PULSE_ENERGY = 100
COLUMNS = ["hWidth", "hHeight", "VzDist", "zDist", "chirp", "b", "pulseEnergy", "intensityMultiplier",
           "hDepth", "hDepthVel", "VxDist", "xDist", "chirpT", "bT", "VzC", "zC", "VxC", "xC"]
INIT_PARAMS = [1.82534513746, 0.00069511761, 0.00027392079, 0.71930270898, 0.00035, 2413.46744543, PULSE_ENERGY,
               1.0, 0.60844837915, 0.00494641054, 0.00493057439, 0.60650040415, 0.00065, 9.83513928219, 0.0, 0.0, 0.0, 0.0]

INIT_PS = pd.DataFrame([INIT_PARAMS], columns=COLUMNS).iloc[0]
# NN Limit is the force (without constant) between two initial_splits/10
#   is valued at 75.032594
NN_LIMIT = pow(INIT_PS["pulseEnergy"]/SPLIT_NUM, 2) / \
    (pow(2*INIT_PS["hHeight"]/SPLIT_NUM, 2) +
     pow(2*INIT_PS["hWidth"]/SPLIT_NUM, 2))/50
