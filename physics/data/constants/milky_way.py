import scipy.constants as spc
import astropy.constants.iau2015 as ac
import astropy.units as u

observedGalacticVelocity = [(229 - 0.2) * 1000, (229 + 0.2) * 1000]
solarObitalVelocity = (observedGalacticVelocity[0] + observedGalacticVelocity[1]) / 2
solarCircleEnclosedMass = 5 * 10**10 * ac.M_sun.value
galacticRadiusObservational = [
    (8.122 - 0.031) * spc.kilo * spc.parsec,
    (8.122 + 0.031) * spc.kilo * spc.parsec,
]
