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

location = EarthLocation.of_site('Roque de los Muchachos') # <--- input location here

# location = EarthLocation(lat=53*deg, lon=1.5*deg, height=131*m) # or manually here


### OBSERVING DATE
observing_date = Time('2025-02-05', location=location) # <--- input observing date here (YYYY-MM-DD)
today = Time(Time.now().iso.split()[0], location=location) # <--- if you want to see tonight

if str(sys.argv[0]) == 'now' or str(sys.argv[0]) == 'today' or str(sys.argv[0]) == 'tonight':
    observing_date = today

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




# SUNRISE AND SUNSET
# sunrise
ax.vlines(sunriseset[0], 0, 90, 
          color='red', 
          linestyle='dashed', 
          alpha=0.75
         ) # line showing sunset

ax.annotate(f'Sunset: {hourText(sunriseset[0])}', 
            (sunriseset[0] + 0.0075*u.d, 71), 
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
            (sunriseset[1] - 0.012*u.d, 71), 
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
            (twlt[0][0] - 0.012*u.d, 71), 
            color='k', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )

ax.annotate(f'Civil: {hourText(twlt[0][1])}', 
            (twlt[0][1] + 0.0075*u.d, 71), 
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
            (twlt[1][0] - 0.012*u.d, 71), 
            color='#909090', 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )

ax.annotate(f'Nauticle: {hourText(twlt[1][1])}', 
            (twlt[1][1] + 0.0075*u.d, 71), 
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
            (twlt[2][0] - 0.012*u.d, 71), 
            color='white', 
            alpha=0.7, 
            rotation=90, 
            fontsize=8, 
            fontweight='bold'
           )

ax.annotate(f'Astronomical: {hourText(twlt[2][1])}', 
            (twlt[2][1] + 0.0075*u.d, 71), 
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
