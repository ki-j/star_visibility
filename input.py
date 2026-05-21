from astropy.coordinates import SkyCoord, EarthLocation, get_body


##############################
# INPUT YOUR PARAMETERS HERE #
##############################

# Welcome to input.py!
# --------------------

# In this file you can input the following information:
#   - Location : either as latitude, longitude, and elevation, or by observatory name (e.g. 'Roque de los Muchachos')
#   - Observing date : either a date in the form 'YYYY-MM-DD' or 'now'
#   - Target(s) : enter name, right ascension, and declination in the form 'XXhXXmXX.XXXs' and 'XXdXXmXX.XXXs', you may enter as many as you like
#   - Telescope : choose a telescope from the list below and set its value to True to plot its shutter limits (or add them yourself)




# LOCATION
# --------

def todeg(h, m, s):
    """ Convert from sexagesimal to degrees. Use for the inputs below. """
    return h+m/60+s/3600

# roque de los muchachos (ORM)
latitude = None # in degrees (float)
longitude = None # in degrees (float)
elevation = None # in metres (float)
observatory = 'Roque de los Muchachos' # e.g. 'Roque de los Muchachos' (string)


# UBO
# observatory = None
# def todeg(h, m, s):
#     return h+m/60+s/3600
# latitude = todeg(52, 23, 14.8)
# longitude = todeg(1, 56, 39.1)
# elevation = 187.7




# OBSERVING DATE
# --------------

# Observing_Date = 'now' # enter a date (YYYY-MM-DD) or the string 'now'
Observing_Date = '2026-05-22' # enter a date (YYYY-MM-DD) or the string 'now'
    

    

# TARGET(S)
# ---------

# test stars
starname = ['Vega', 'V1315 Aquilae', 'LAMOST J035913.61+405035.0'] # enter star name(s) into this list (e.g. 'Vega')
RA = ["18h36m56.33635s", "19h13m54.5308677240s", "03h59m13.626515992s"] # enter right ascension value(s) here in the form 'XXhXXmXX.XXXs'
DEC = ["+38d47m01.2802s", "+12d18m03.239745228s", "+40d50m35.095271748s"] # enter declination value(s) here in the form 'XXdXXmXX.XXXs'

# Issy's GRB (Mar 2026)
# starname = ["Issy's GRB"]
# RA = ["14h14m00s"]
# DEC = ["78d42m00s"]




# TELESCOPE
# ---------

INT = False
TNG = True


###############################################################