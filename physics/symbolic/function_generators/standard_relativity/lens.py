from physics.symbolic.variables import G, c, M, b

def gravitationalLensDeflection(b = b):
    return 4 * G * M / c**2 / b
