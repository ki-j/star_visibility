# imports
import sys
import numpy as np
from matplotlib import pyplot as plt

from astropy.time import Time
from astropy import units as u
from astropy import constants as const
from astropy.visualization import time_support, quantity_support
quantity_support()
time_support()
from astropy.coordinates import SkyCoord, EarthLocation, get_body

from helper_functions import whatsGoingOnTonight, hourText, find_moondist

### INPUTS #################################################

# todo : put location info through as an argument for terminal use

### LOCATION

# location = EarthLocation.of_site('Roque de los Muchachos') # <--- input location here

location = EarthLocation(lat=52*u.deg, lon=-1.8*u.deg, height=140*u.m) # <-- birmingham

# location = EarthLocation(lat=53*deg, lon=1.5*deg, height=131*m) # or manually here


### OBSERVING DATE
obsdate = Time('2025-02-05', location=location) # <--- input observing date here (YYYY-MM-DD) (if not looking for current observing conditions)
    
now = Time(Time.now(), location=location, format='iso')
today = Time(now.iso.split()[0], location=location) # <--- if you want to see tonight

if str(sys.argv[1]) == 'now' or str(sys.argv[1]) == 'today' or str(sys.argv[1]) == 'tonight':
    if float(now.datetime.hour) > 12.: # becaue of how the day is defined. Sunset is the date - 1 day, easier just to add a day now so I did
        observing_date = today  + 1*u.day # vis today
    else:
        observing_date = today
    plotnow = True
else:
    observing_date = obsdate + 1*u.day # vis on given date
    plotnow = False
    

### TARGET(S)
vega = SkyCoord(ra="18h36m56.33635s", dec="+38d47m01.2802s", frame="icrs")
V1315aql = SkyCoord(ra="19h13m54.5308677240s", dec="+12d18m03.239745228s", frame="icrs")
LAMOST0359 = SkyCoord(ra="03h59m13.626515992s", dec="+40d50m35.095271748s", frame="icrs")


# ADD TO LISTS HERE ALSO
targets = [vega, V1315aql, LAMOST0359]
starname = ['vega', 'V1315aql', 'LAMOST0359']

# SUN AND MOON
sun = get_body('sun', observing_date)
moon = get_body('moon', observing_date)

###############################################################


### UNITS #####################################################

deg = u.degree
m = u.meter

###############################################################


### RUN STARALT ###############################################

# calculate all stuff for moon
night_moon, altitudes_moon, azimuths_moon, sunriseset_moon, twlt_moon, trans_moon, midnight_moon = whatsGoingOnTonight(location, observing_date, moon)


# calculate all stuff for target(s)
target_list = []
separations = [] # separation from moon

for target in targets:
    
    night, altitudes, azimuths, sunriseset, twlt, trans, midnight = whatsGoingOnTonight(location, observing_date, target)
    target_list.append([night, altitudes])
    
    # separation from moon
    seps_target = []
    
    for i in range(len(altitudes)):
        mn = [altitudes_moon[i], azimuths_moon[i]]
        targ = [altitudes[i], azimuths[i]]
        seps_target.append(find_moondist(mn, targ))
        
    separations.append(seps_target)
    
        
# night, altitudes, azimuths, sunriseset, twlt, trans, midnight = whatsGoingOnTonight(location, observing_date, targets[0])
night = target_list[0][0] # array of times for the night


# create the plot
fig, ax = plt.subplots(figsize=[20, 10])


textHeight = 65 # the bottom of the annotation text (in degrees of altitude)

# SUNRISE AND SUNSET
# sunrise
ax.vlines(sunriseset[0], 0, 90, 
          color='red', 
          linestyle='dashed', 
          alpha=0.75
         ) # line showing sunset

ax.annotate(f'Sunset: {hourText(sunriseset[0])}', 
            (sunriseset[0] + 0.0075*u.d, textHeight), 
            color='red', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )


# sunset
ax.vlines(sunriseset[1], 0, 90, 
          color='red', 
          linestyle='dashed', 
          alpha=0.75
         ) # line showing sunrise

ax.annotate(f'Sunrise: {hourText(sunriseset[1])}', 
            (sunriseset[1] - 0.012*u.d, textHeight), 
            color='red', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )

# fill
ax.fill_between(np.linspace(sunriseset[0], sunriseset[1]), 0, 90, 
                color='#000510', 
                alpha=0.4
               ) # darken




# TWILIGHT
# civil
ax.fill_between(np.linspace(twlt[0][0], twlt[0][1]), 0, 90, 
                color='#000510', 
                alpha=0.5
               )

ax.annotate(f'Civil: {hourText(twlt[0][0])}', 
            (twlt[0][0] - 0.012*u.d, textHeight), 
            color='k', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )

ax.annotate(f'Civil: {hourText(twlt[0][1])}', 
            (twlt[0][1] + 0.0075*u.d, textHeight), 
            color='k', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )


# nauticle
ax.fill_between(np.linspace(twlt[1][0], twlt[1][1]), 0, 90, 
                color='#000510', 
                alpha=0.5
               )

ax.annotate(f'Nauticle: {hourText(twlt[1][0])}', 
            (twlt[1][0] - 0.012*u.d, textHeight), 
            color='#909090', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )

ax.annotate(f'Nauticle: {hourText(twlt[1][1])}', 
            (twlt[1][1] + 0.0075*u.d, textHeight), 
            color='#909090', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )


# astronomical
ax.fill_between(np.linspace(twlt[2][0], twlt[2][1]), 0, 90, 
                color='#000510', 
                alpha=0.5
               )

ax.annotate(f'Astronomical: {hourText(twlt[2][0])}', 
            (twlt[2][0] - 0.012*u.d, textHeight), 
            color='white', 
            alpha=0.7, 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )

ax.annotate(f'Astronomical: {hourText(twlt[2][1])}', 
            (twlt[2][1] + 0.0075*u.d, textHeight), 
            color='white', 
            alpha=0.7, 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )


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

# vertical 'now' line
if plotnow:
    ax.vlines(now, 0, 90,
             color='cyan',
             linestyle='dashed',
             alpha=1.
                 ) # plot vertical 'now' line

    ax.annotate(f'{hourText(now)}', 
            (now - 0.012*u.d, textHeight), 
            color='cyan', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           ) # annotate now line


# LIMITS, LABELS & TITLE
ax.set_title(f'{night[0].to_string().split()[0]} $|$ lat={location.lat:.2f} long={location.lon:.2f} elev={location.height:.0f}')

labels = []
tickvals = []
for i in range(len(night)-1):
    if night[i].datetime.hour < night[i+1].datetime.hour or night[i].datetime.hour > night[i+1].datetime.hour:
        labels.append(night[i].datetime.hour)
        tickvals.append(night[i]-1*u.hour)

# import pdb; pdb.set_trace()

ax.set_xticks(ticks=Time(tickvals))
ax.set_xticklabels(labels)

ax.set_ylim(0, 90), ax.set_ylabel('Altitude (degrees)')

ax.set_xlim(sunriseset[0]-0.75*u.hour, sunriseset[-1]+0.75*u.hour)

plt.legend()

fig.savefig('STARALT.pdf', bbox_inches='tight')

plt.show()
