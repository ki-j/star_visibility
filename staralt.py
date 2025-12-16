# imports
import numpy as np
from matplotlib import pyplot as plt

# import argparse

from astropy.time import Time
from astropy import units as u
from astropy import constants as const
from astropy.visualization import time_support, quantity_support
quantity_support()
time_support()
from astropy.coordinates import SkyCoord, EarthLocation, get_body

# locations



# parser
# parser = argparse.ArgumentParser()
# parser.add_argument("location")
# parser.add_argument("obsDate")
# parser.add_argument("target")
# args = parser.parse_args()

# inputs from command line
# location
# if type(args.location) == str: # if inputting a location valid with astropy.coordinates
#     location = EarthLocation.of_site(str(args.location))
# if type(args.location) == list: # if inputting the latitude, longitude, and elevation as a list
#     location = EarthLocation(lat=args.location[0]*deg, lon=args.location[1]*deg, height=args.location[2]*m)
# else: # if inputting the full EarthLocation
#     location = args.location
# if type(location) != 'astropy.coordinates.earth.EarthLocation' and type(location) != str and type(location) != list:
#     raise TypeError('Location has been incorrectly formatted. Please include either a valid location, a list in the form [lat, long, height], or an EarthLocation() object.')

# obsDate
# if type(args.obsDate) != str:
#     raise TypeError('Observing date should be a string in the form YYYY-MM-DD')
# else:
#     obsDate = Time(args.obsDate, location=location)

# target
# if type(args.target) != list:
#     raise TypeError('Please write the target as a list in the form [RAhRAmRAs, +DECdDECmDECs]')
# else:
#     target = SkyCoord(ra=str(args.target[0]), dec=str(args.target[1]), frame="icrs")

# units
deg = u.degree
m = u.meter
    
# calculate transit time and altitude for a given target
def transits(location, target, time):
    """
    determine the transit time, altitude and middle of the night of a target given a location
    """
    # transit altitude
    transitAltitude = (90*deg - location.lat) + target.dec
    # RA of the Sun on the observing date
    raSun = get_body('sun', time).ra.hour *u.hour
    # middle of the nigth (UTC)
    timeMidnight = -location.lon.hour * u.hour + time
    # optimal RA of the Sun
    def optRAsun(target):
        raSunOptimal = target.ra.hour*u.hour + 12*u.hour
        if raSunOptimal > 24*u.hour:
            raSunOptimal = raSunOptimal - 24*u.hour
        return raSunOptimal
    raSunOptimal = optRAsun(target)
    # difference between optimal and observing date RA_sun
    dRA = raSunOptimal - raSun
    # find the transit time
    transitTime = timeMidnight + dRA
    return transitAltitude.deg, transitTime, timeMidnight

# determine sunrise and sunset times
def riseset(location, time):
    # get sun
    sun = get_body('sun', time)
    # find hour angle of sunset
    numerator = np.sin(-0.5 *np.pi/180) - np.sin(location.lat.rad) * np.sin(sun.dec.rad)
    denominator = np.cos(location.lat.rad) * np.cos(sun.dec.rad)
    haSun = np.arccos(numerator/denominator) *180/(15*np.pi)
    # find sunset
    sunset = 12*u.hour + haSun*u.hour + time - 1*u.day
    # difference between sunset and middle of night
    midnight = -location.lon.hour * u.hour + time
    time2midnight = abs(sunset - midnight).value * 24 * u.hour
    # sunrise
    sunrise = midnight + time2midnight
    return sunset, sunrise

# calculate the twilight times
def twilight(location, time, angle):
    """
    calculate the times of twilight for a given night & location and angle
    """
    # sun
    sun = get_body('sun', time)
    # find hour angle of twilight
    numerator = np.sin(-angle *np.pi/180) - np.sin(location.lat.rad) * np.sin(sun.dec.rad)
    denominator = np.cos(location.lat.rad) * np.cos(sun.dec.rad)
    haSun = np.arccos(numerator/denominator) *180/(15*np.pi)
    # find end of twilight
    twilight_e = 12*u.hour + haSun*u.hour + time - 1*u.day
    # difference between twilight and middle of night
    midnight = -location.lon.hour * u.hour + time
    time2midnight = abs(twilight_e - midnight).value * 24 * u.hour
    # start of twilight
    twilight_s = midnight + time2midnight
    return twilight_s, twilight_e

def twilight_zones(location, time):
    """
    calculate the start and end times for civil, nauticle and astronomical twillight
    """
    # civil twilight
    civil = twilight(location, time, 6) # sun 6 degrees below horizon
    # nauticle twilight
    nauticle = twilight(location, time, 12) # sun 12 degrees below horizon
    # civil twilight
    astronomical = twilight(location, time, 18) # sun 18 degrees below horizon
    return civil, nauticle, astronomical

# calculate hour angle
def find_ha(location, time, target, unit='degree'):
    """
    calculate the hour angle of a given target for some location and time
    """
    # find LST
    LST = time.sidereal_time('apparent')
    # calculate HA
    HA = LST.deg - target.ra.deg
    
    if unit == 'hour':
        HA = (HA / 15) *u.hour # return HA in hours
    else:
        HA = HA *u.deg # return HA in degrees (default)
    return HA

# calculate altitude
def find_alt(location, target, ha):
    lat = location.lat.rad # latitude in radians
    dec = target.dec.rad # declination in radians
    if type(ha) != list:
        sinalt = np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(ha.to('radian'))
        alt = np.arcsin(sinalt)
        alts = alt.to(deg)
    else:
        alts = []
        for h in ha:
            sinalt = np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(h.to('radian'))
            alt = np.arcsin(sinalt)
            alts.append(alt.to(deg))
    return alts

# determine the full span of the night, depending on sunrise and set times
def find_night(location, obs_date):
    """
    create an array of times for a given date to define the night (+/- hours before/after midnight)
    """
    # find sunrise and sunset
    srise, sset = riseset(location, obs_date)
    # return an array of times per minute
    return np.linspace(srise - 1*u.hour, sset + 1*u.hour, 24*60)

# define the target's trajectory over the night
def whatsGoingOnTonight(location, obs_date, target):
    """
    determine a target's altitude, transit times across a night
    """
    # night
    night = find_night(location, obs_date) # define the observing night
    sset, srise = riseset(location, obs_date) # define the sunrise and sunset for the goven night
    civil, nauticle, astronomical = twilight_zones(location, obs_date) # define civil, nauticle, astronomical twilight
    
    # target
    altitudes = [] # array of altitudes as a function of time
    for minute in night:
        ha = find_ha(location, minute, target, unit='degree')
        alt = find_alt(location, target, ha)
        altitudes.append(alt.value)
    transAlt, transTime, mid = transits(location, target, obs_date) # transit altitude and time
    
    return night, altitudes, [sset, srise], [civil, nauticle, astronomical], [transAlt, transTime], mid

def hourText(time):
    """
    take the Time() object and return just hh:mm UTC
    """
    t = str(time).split()
    hr = t[1].split(":")[0] + ':' + t[1].split(":")[1] + ' UTC'
    return hr

# full plotting function, including Moon trajectory
# turn it into a function
def plotAlt(location, obsDate, target, fs=[9,10], starname='Star'):
    """
    take a location, observing date and target and plot altitude over the night
    """

    # objects
    sun = get_body('sun', obsDate)
    moon = get_body('moon', obsDate)
    
    # calculate all stuff for target
    night, altitudes, sunriseset, twlt, trans, midnight = whatsGoingOnTonight(location, obsDate, target)
    # calculate all stuff for moon
    night_moon, altitudes_moon, sunriseset_moon, twlt_moon, trans_moon, midnight_moon = whatsGoingOnTonight(location, obsDate, moon)
    
    # create the plot
    fig, ax = plt.subplots(figsize=fs)
    
    # SUNRISE AND SUNSET
    # sunrise
    ax.vlines(sunriseset[0], 0, 90, 
              color='red', linestyle='dashed', alpha=0.75) # line showing sunset
    ax.annotate(f'Sunset: {hourText(sunriseset[0])}', (sunriseset[0]+0.0075, 71), 
                color='red', rotation=90, fontsize=8, fontweight='bold')
    # sunset
    ax.vlines(sunriseset[1], 0, 90, 
              color='red', linestyle='dashed', alpha=0.75) # line showing sunrise
    ax.annotate(f'Sunrise: {hourText(sunriseset[1])}', (sunriseset[1]-0.012, 71), 
                color='red', rotation=90, fontsize=8, fontweight='bold')
    # fill
    ax.fill_between(np.linspace(sunriseset[0], sunriseset[1]), 0, 90, 
                    color='k', alpha=0.4) # darken
    
    # TWILIGHT
    # civil
    ax.fill_between(np.linspace(twlt[0][0], twlt[0][1]), 0, 90, color='k', alpha=0.5)
    ax.annotate(f'Civil: {hourText(twlt[0][0])}', (twlt[0][0]-0.012, 71), 
                color='k', rotation=90, fontsize=8, fontweight='bold')
    ax.annotate(f'Civil: {hourText(twlt[0][1])}', (twlt[0][1]+0.0075, 71), 
                color='k', rotation=90, fontsize=8, fontweight='bold')
    # nauticle
    ax.fill_between(np.linspace(twlt[1][0], twlt[1][1]), 0, 90, color='k', alpha=0.5)
    ax.annotate(f'Nauticle: {hourText(twlt[1][0])}', (twlt[1][0]-0.012, 71), 
                color='#909090', rotation=90, fontsize=8, fontweight='bold')
    ax.annotate(f'Nauticle: {hourText(twlt[1][1])}', (twlt[1][1]+0.0075, 71), 
                color='#909090', rotation=90, fontsize=8, fontweight='bold')
    # astronomical
    ax.fill_between(np.linspace(twlt[2][0], twlt[2][1]), 0, 90, color='k', alpha=0.5)
    ax.annotate(f'Astronomical: {hourText(twlt[2][0])}', (twlt[2][0]-0.012, 71), 
                color='white', alpha=0.7, rotation=90, fontsize=8, fontweight='bold')
    ax.annotate(f'Astronomical: {hourText(twlt[2][1])}', (twlt[2][1]+0.0075, 71), 
                color='white', alpha=0.7, rotation=90, fontsize=8, fontweight='bold')
    
    # ALTITUDE
    ax.plot(night, altitudes, color='orange', alpha=0.75, label=starname)
    ax.plot(night_moon, altitudes_moon, color='grey', linestyle='dashed', alpha=0.75, label='Moon')

    # LIMITS, LABELS & TITLE
    ax.set_title(starname)
    labels = []
    hour = night[0].datetime.hour
    end = night[-1].datetime.hour + 1
    while hour != end:
        if hour > 23:
            hour = hour - 24
        labels.append(hour)
        hour += 1
        if hour == end:
            break

    ax.set_xticks(np.linspace(night[0], night[-1], len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 90), ax.set_ylabel('Altitude (degrees)')
    ax.set_xlim(night[0], night[-1])
    plt.legend()
    plt.show()

    print(f"Transit time: {hourText(trans[1])}")
    print(f"Transit altitude {trans[0]:.1f}")
    
    return trans


# plotAlt(location, obsDate, target)