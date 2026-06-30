import numpy as np

def superFactorial(n, step = 1):
    arr = np.arange(0, n, step).cumsum()
    return arr[arr.__len__() - 1]
