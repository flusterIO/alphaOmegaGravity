import numpy as np


def eta():
    return np.array([[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


def metricTensor(x0=-1, x1=1, x2=1, x3=1):
    return np.array([[x0, 0, 0, 0], [0, x1, 0, 0], [0, 0, x2, 0], [0, 0, 0, x3]])

# def metricTensorPolar():

