# imports
import numpy as np
from matplotlib import pyplot as plt

from astropy.time import Time
from astropy import units as u
from astropy import constants as const
from astropy.visualization import time_support, quantity_support
quantity_support()
time_support()
from astropy.coordinates import SkyCoord, EarthLocation, get_body

### LOCATIONS #################################################

# todo : put location info through as an argument for terminal use
location = EarthLocation.of_site('Roque de los Muchachos') # <--- input location here
observing_date = Time('2025-02-05', location=location) # <--- input observing date here
target = SkyCoord(ra='05h55m10.30536s', dec='+07d24m25.4304s', frame='icrs') # <--- input target here (e.g. betelgeuse)
sun = get_body('sun', observing_date)
moon = get_body('moon', observing_date)

###############################################################


### UNITS #####################################################

deg = u.degree
m = u.meter

###############################################################


### TRANSITS ##################################################

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

###############################################################


### SUNRISE/SET & TWILIGHT ####################################
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

###############################################################


### HOUR ANGLE ################################################

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

###############################################################


### ALTITUDE ##################################################

def find_alt(location, target, ha):
    """ find the altitude of an objects for a given location and altitude """
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

###############################################################


### AZIMUTH ###################################################

def find_az(location, target, altitude):
    """ calculate the azimuth of a target at a given time given the altitude (and therefore time) from spherical trig """
    lat = location.lat.rad # latitude in radians
    dec = target.dec.rad # declination in radians
    alt = altitude # altitude in degrees
    # from spherical trig
    num = np.sin(dec) - np.sin(lat)*np.sin(alt.to('radian'))
    den = np.cos(lat)*np.cos(alt.to('radian'))
    cosaz = num / den
    az = np.arccos(cosaz)
    if type(az) == float:
        az = az * u.radian
    return az.to(deg)

###############################################################


### NIGHT  ####################################################

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
    azimuths = []
    for alt in altitudes:
        alt = alt * u.deg
        azimuths.append(find_az(location, target, alt).value)
    
    return night, altitudes, azimuths, [sset, srise], [civil, nauticle, astronomical], [transAlt, transTime], mid

###############################################################


### MOON SEPARATION ###########################################

# calculate Moon distance
def find_moondist(target1, target2, ang='degrees'):
    """ calculate the angular separation between two objects in the night sky (or the daytime sky, I won't judge) 

    params
    ------
    target1, 2 (list [float, float]) : altaz coordinates of targets in the format [alt, az], in degrees preferably (if not change ang to 'silly bastard')

    output
    ------
    separation (float) : angular separation of the two targets in degrees
    """
    # params
    alt1, az1 = target1 # in degrees remember!
    alt2, az2 = target2
    # convert to equation form
    z1 = 90 - alt1
    z2 = 90 - alt2
    gamma = np.abs(az1 - az2)
    if gamma > 180: # don't go round the whole circle like a dummy!
        gamma -= 180
    # convert to radians for calculations
    z1 *= np.pi/180
    z2 *= np.pi/180
    gamma *= np.pi/180
    # cosine rule
    cosAB = np.cos(z1)*np.cos(z2) - np.sin(z1)*np.sin(z2)*np.cos(gamma)
    AB = np.arccos(cosAB)
    # back to degrees
    separation = AB * 180/np.pi
    return separation

###############################################################


### PLOTTING ##################################################

def hourText(time):
    """
    take the Time() object and return just hh:mm UTC
    """
    t = str(time).split()
    hr = t[1].split(":")[0] + ':' + t[1].split(":")[1] + ' UTC'
    return hr

# full plotting function
# turn it into a function
def plotAlt(location, obsDate, targets, fs=[20,10], starname=['Star']):
    """
    take a location, observing date and target(s) and plot altitude over the night
    """
    # calculate all stuff for moon
    night_moon, altitudes_moon, azimuths_moon, sunriseset_moon, twlt_moon, trans_moon, midnight_moon = whatsGoingOnTonight(location, obsDate, moon)
    
    # calculate all stuff for target
    target_list = []
    separations = [] # separation from moon
    for target in targets:
        night, altitudes, azimuths, sunriseset, twlt, trans, midnight = whatsGoingOnTonight(location, obsDate, target)
        target_list.append([night, altitudes])
        # separation from moon
        seps_target = []
        for i in range(len(altitudes)):
            mn = [altitudes_moon[i], azimuths_moon[i]]
            targ = [altitudes[i], azimuths[i]]
            seps_target.append(find_moondist(mn, targ))
        separations.append(seps_target)
            
    night, altitudes, azimuths, sunriseset, twlt, trans, midnight = whatsGoingOnTonight(location, obsDate, targets[0])
    # calculate all stuff for moon
    night_moon, altitudes_moon, azimuths_moon, sunriseset_moon, twlt_moon, trans_moon, midnight_moon = whatsGoingOnTonight(location, obsDate, moon)
    
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
                    color='#000510', alpha=0.4) # darken
    
    # TWILIGHT
    # civil
    ax.fill_between(np.linspace(twlt[0][0], twlt[0][1]), 0, 90, color='#000510', alpha=0.5)
    ax.annotate(f'Civil: {hourText(twlt[0][0])}', (twlt[0][0]-0.012, 71), 
                color='k', rotation=90, fontsize=8, fontweight='bold')
    ax.annotate(f'Civil: {hourText(twlt[0][1])}', (twlt[0][1]+0.0075, 71), 
                color='k', rotation=90, fontsize=8, fontweight='bold')
    # nauticle
    ax.fill_between(np.linspace(twlt[1][0], twlt[1][1]), 0, 90, color='#000510', alpha=0.5)
    ax.annotate(f'Nauticle: {hourText(twlt[1][0])}', (twlt[1][0]-0.012, 71), 
                color='#909090', rotation=90, fontsize=8, fontweight='bold')
    ax.annotate(f'Nauticle: {hourText(twlt[1][1])}', (twlt[1][1]+0.0075, 71), 
                color='#909090', rotation=90, fontsize=8, fontweight='bold')
    # astronomical
    ax.fill_between(np.linspace(twlt[2][0], twlt[2][1]), 0, 90, color='#000510', alpha=0.5)
    ax.annotate(f'Astronomical: {hourText(twlt[2][0])}', (twlt[2][0]-0.012, 71), 
                color='white', alpha=0.7, rotation=90, fontsize=8, fontweight='bold')
    ax.annotate(f'Astronomical: {hourText(twlt[2][1])}', (twlt[2][1]+0.0075, 71), 
                color='white', alpha=0.7, rotation=90, fontsize=8, fontweight='bold')


    # ALTITUDE & MOON SEPARATION
    colours = ['#005F73', '#0A9396', '#94D2BD', '#E9D8A6', '#EE9B00', '#CA6702', '#BB3E03', '#AE2012'] # ignore this if you want lol
    # colours = ['#D8F3DC', '#B7E4C7', '#95D5B2', '#74C69D', '#52B788', '#40916C', '#2D6A4F', '#1B4332'] # green palette
    for i in range(len(target_list)):
        times = target_list[i][0]
        alts = target_list[i][1]
        ax.plot(times, alts, alpha=0.75, color=colours[-i], label=starname[i])
        for k in range(15): # annotate the moon separations
            index = int(k/15 * len(times)) # split times into 15 equal points
            ax.annotate(f"{separations[i][index]:.0f}", (times[index], alts[index]), fontsize=10, color=colours[-i])
    ax.plot(night_moon, altitudes_moon, color='grey', linestyle='dashed', alpha=0.75)
    
    
    # LIMITS, LABELS & TITLE
    ax.set_title(f'{night[0].to_string().split()[0]} $|$ lat={location.lat:.2f} long={location.lon:.2f} elev={location.height:.0f}')
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
    fig.savefig('STARALT.pdf', bbox_inches='tight')
    plt.show()

    # print(f"Transit time: {hourText(trans[1])}")
    # print(f"Transit altitude {trans[0]:.1f}")
    
    return trans

###############################################################

# plotAlt(location, obsDate, target)