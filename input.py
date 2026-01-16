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
##########

latitude = None # in degrees (float)
longitude = None # in degrees (float)
elevation = None # in metres (float)
observatory = 'Roque de los Muchachos' # e.g. 'Roque de los Muchachos' (string)



# OBSERVING DATE
################

Observing_Date = '2000-06-09' # enter a date (YYYY-MM-DD) or the string 'now'
    
    

# TARGET(S)
###########

starname = ['Vega', 'V1315 Aquilae', 'LAMOST J035913.61+405035.0'] # enter star name(s) into this list (e.g. 'Vega')
RA = ["18h36m56.33635s", "19h13m54.5308677240s", "03h59m13.626515992s"] # enter right ascension value(s) here in the form 'XXhXXmXX.XXXs'
DEC = ["+38d47m01.2802s", "+12d18m03.239745228s", "+40d50m35.095271748s"] # enter declination value(s) here in the form 'XXdXXmXX.XXXs'


# TELESCOPE
###########

INT = True


###############################################################
