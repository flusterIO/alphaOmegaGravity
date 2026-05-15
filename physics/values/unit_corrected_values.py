from astropy.constants.constant import Quantity
import astropy.units as u
import astropy.constants.iau2015 as ac
import scipy.constants as spc


unitCorrected_G: Quantity = spc.G * u.m**3 / u.kg / u.s

unitCorrected_g: Quantity = unitCorrected_G * ac.M_earth / ac.R_earth**2

unitCorrected_alpha: Quantity = spc.alpha * u.m
