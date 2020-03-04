import numpy as np


def load(file):
    return np.load(file, allow_pickle=True).item()
