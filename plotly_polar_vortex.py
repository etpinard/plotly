# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# How cold has this winter been?

# <markdowncell>

# Let's check using IPython and plot.ly!
# 
# I will plot a world map of the average daily temperature anomalies (i.e. deviation from climatology) for December 2013 and Februray 2014.
# 
# This notebook will use:
# 
# * plot.ly's 'heatmap' plot type for colored contour maps compatible with plot.ly [see https://plot.ly/api/python/docs/heatmaps ]
# * NetCDF data imports which are very common in Atmospheric/Oceanic sciences [ref. http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.netcdf.netcdf_file.html ]
# * Basemap module for coastlines and countries data [ref. http://matplotlib.org/basemap/ ]
# 
# 
# Note that, I would have loved to use Basemap also to project the data to a curvilinear grid (e.g. the Robinson projection for aesthestics). <br>
# However, plot.ly's 'heatmap' plot type currently only accepts 1-dimensional 'x' and 'y' coordinate vectors. <br>
# So, my plot will use the simplest of all cylindrical projections (called the equirectangular projection) where all meridians have constant spacing.
# 
# I divided the tasks as:
# 
# - Setup plot.ly
# - Import the temperature data
# - Get the coastline and country boundary data
# - Setup plot layout and color map
# - Plot!

# <headingcell level=2>

# 1) Setup plot.ly

# <markdowncell>

# Import plot.ly with python api [see https://plot.ly/api/python/getting-started ], 
# sign in, create shortcut to plot.ly object and check version:

# <codecell>

import plotly

username='etpinard'       # my username
api_key = 'a35m7g6el5'    # my api key

py = plotly.plotly(username,api_key)  # shortcut to plot.ly object
plotly.__version__                    # check ploty.ly version

# <markdowncell>

# Define a function to embed plot.ly graph in IPython notebook [ref. http://nbviewer.ipython.org/github/plotly/IPython-plotly/blob/master/Plotly%20Quickstart.ipynb ]:

# <codecell>

from IPython.display import HTML

def show_plot(url, width=1000, height=500):
    s = '<iframe height="%s" id="igraph" scrolling="no" seamless="seamless" src="%s" width="%s"></iframe>' %\
    (height+50, "/".join(map(str,[url, width, height])), width+50)
    return HTML(s)

# <markdowncell>

# Import all other packages needed for executing this notebook:

# <codecell>

import numpy as np                         # numpy ... of course
from scipy.io import netcdf                # scipy NetCDF i/o module
from mpl_toolkits.basemap import Basemap   # Basemap model of the matplotlib toolkit

# <headingcell level=2>

# 2) Import temperature data

# <markdowncell>

# I got data on this http://www.esrl.noaa.gov/psd/data/composites/day/ .
# 
# This website does not allow to *code* your output demand and/or use `wget` to download the data. <br>
# However, the data I used can be found in a only a few clicks:
# 
# - Select `Air Temperature` in `Varaibles`
# - Select `Surface` in `Analysis level?`
# - Select `Dec` | `1` and `Jan` | `31` 
# - Enter `2014` in the `Enter Year of last day of range` prompt
# - Select `Anomaly` in `Plot type?`
# - Select `All` in `Region of globe`
# - Click on `Create Plot` and a temporary NetCDf file can be downloaded on that new.

# <markdowncell>

# The data represents the Dec2013-Jan2014 average daily surface air temperature anomaly (in deg. C) with respect to 1981-2010 climatology.
# 
# That is, our temperature variable of interest $T$ is
# 
# $$ T = \frac{1}{62 \; \text{days}} \sum_{i=1}^{62} \bigg( T_i - \bar{\hskip2pt T}_i \bigg) $$
# 
# where $\bar{\hskip2pt T}_i$ is daily surface air temepature 1981-2010 climatological mean
# and $i$ is the day index starting from 1 December 2013.
# 
# In the NetCDF file, $T$ is refered to as `air`. <br>
# The NetCDF file also includes a longitude (`lon`), latitude (`lat`) and time singleton (`time`). <br>
# Note also that this dataset is of quite low resolution 2.5 deg x 2.5 deg (about 300 km) resolution.
# 
# Now, let's extract these arrays [inspired by http://earthpy.org/interpolation_between_grids_with_basemap.html ]:

# <codecell>

f_tas_path = 'compday.qQbv_4_mY1.nc'         # file must be in working directory
f_tas = netcdf.netcdf_file(f_tas_path, 'r')  # will most likely have a different name if you downloaded it yourself

lon = f_tas.variables['lon'][::]             # copy as an array
lat = f_tas.variables['lat'][::-1]           # invert the latitude vector, now South to North
air = f_tas.variables['air'][0,::-1,:]       # squeeze out the time dimension, invert latitude index

# shift 'lon' from [0,360] to [-180,180] for better-looking plots
lon_tmp = lon.copy()
for n, l in enumerate(lon_tmp):
    if l >= 180:
        lon_tmp[n] = lon_tmp[n] - 360. 
lon = lon_tmp                         # gives [0,180]U[-180,-2.5]
mid_lon = lon.shape[0]/2              # index of second of the -180 lon
lon_east = lon[0:mid_lon]             # [0,180]
lon_west = lon[mid_lon:]              # [-180,0]
lon = np.hstack((lon_west, lon_east)) # stack the 2 halves

# correspondingly shift the 'air' array
air_east = air[:,0:mid_lon]
air_west = air[:,mid_lon:]
air = np.hstack((air_west, air_east))

f_tas.close()   # close NetCDF file

#print lon
#print lat

# <headingcell level=2>

# 3) Get coastline and country boundary data

# <markdowncell>

# The Basemap module includes data for drawing coastlines and country boundaries onto world maps. <br>
# Adding coastlines and/or country boundaries on a regular Python (i.e, non-plot.ly) figure can easily be done
# with the `.drawcoaslines()` or `.drawcountries()` Basemap methods. 
# 
# Here, I retrive the Basemap plotting data (or polygons) and convert them to longitude/latitude arrays compatible with plot.ly <br>
# [inspired by http://stackoverflow.com/questions/14280312/world-map-without-rivers-with-matplotlib-basemap ]
# 
# In other words, the goal is to plot each *continuous* coastline and country boundraies as 1 plot.ly *trace*.
# 
# Note that there are other ways to this, such as import shapefile file off the web [e.g. http://www.naturalearthdata.com/downloads/110m-cultural-vectors/ ].

# <codecell>

m = Basemap()     # shortcut to Basemap object (no need to specify projection type)

# Function assigning plot.ly plotting options for each coastline/country trace
def MakePlotData(lon_cc,lat_cc):
    ''' lon_cc: lon. of coastline/country, lat_cc: lat of coastline/country '''
    return {
            'x': lon_cc,
            'y': lat_cc,  
            'mode': 'lines',
            'line': {'color': 'rgb(0,0,0)'},
            'name': ''
           }

# Functions converting coastline/country plot polygons to lon/lat arrays
def polygons_to_lonlat(N_poly,poly_paths):
    ''' N_poly: number of polygon to convert, poly_paths: paths to polygons '''
    data_tmp = []  # init. plotting list 
    for i_poly in range(N_poly):
        poly_path = poly_paths[i_poly]
        
        # get the Basemap coordinates of each segment
        coords_cc = np.array([(vertex[0],vertex[1]) 
                              for (vertex,code) in poly_path.iter_segments(simplify=False)])
        
        # convert coordinates to lon/lat by 'inverting' the Basemap projection
        lon_cc, lat_cc = m(coords_cc[:,0],coords_cc[:,1], inverse=True)
        
        # add plot.ly plotting options
        data_tmp.append(MakePlotData(lon_cc,lat_cc))
        
    return data_tmp
        
# 1) 'draw' coastlines using the Basemap built-in methods
m_draw = m.drawcoastlines()
poly_paths = m_draw.get_paths()  # get coastline polygon paths
N_poly = 91                      # use of the 91 biggest coastlines (i.e. don't show rivers)
data_cc1 = polygons_to_lonlat(N_poly,poly_paths)

# 2) Similarly 'draw' countries boundaries
m_draw = m.drawcountries()
poly_paths = m_draw.get_paths() # get country boundaries polygon paths
N_poly = len(poly_paths)     # use all countries boundaries
data_cc2 = polygons_to_lonlat(N_poly,poly_paths)

# Concatenate coastlines and country boundaries list
data_cc = data_cc1+data_cc2

#print data_cc
#print len(data_cc), len(data_cc1), len(data_cc2)

# <headingcell level=2>

# 4) Setup plot layout and color map

# <markdowncell>

# Layout options:

# <codecell>

layout_style = {
                'autosize': False,
                'width': 1000,
                'height': 500,
                'margin': {'t':100, 'b':50, 'l':50, 'r':0, 'pad':4},
                "plot_bgcolor": "rgb(255,255,255)",                     # white plot bg
                "paper_bgcolor": "rgb(222,222,222)",                    # gray frame bg
                'showlegend': False,
                'hovermode': 'closest',
                'title': "Mean daily surface air temperature anomalies [in deg. C] <br>\
                          from 2013-12-01 to 2014-01-31 with respect to 1981-2010 climatology",   
                'titlefont': {'color':'rgb(0,0,0)', 'size':18},
                'xaxis': {
                          #'title': 'longitude',
                          'range': [lon[0],lon[-1]],       # show the whole globe
                          'showgrid': False,
                          'zeroline': None,
                          'ticks': None,
                          'showticklabels': False
                         },
                'yaxis': {
                          #'title': 'latitude',
                          'range': [lat[0],lat[-1]],        # show the whole globe
                          'showgrid': True,
                          'zeroline': False,
                          'ticks': None,
                          'showticklabels': False
                         }
                }

# <markdowncell>

# By default, plot.ly linearly maps [`min(z)`, `max(z)`] to [0,1] to build its color map. <br>
# Instead, let's build a symmetric color map (and color bar) around zero, 
# for a more elegant plot.
# 
# So, first find the smaller possible symmetric range with integer bounds 
# that includes all `z` values. <br>
# Then, using the default 'scl' map (coded below as `default_scl_map()`), 
# map this range for our plot.ly call.
# 
# To put emphasis on the temperature anomalies, 
# let's use a variant of ColorBrewer2.0's divergent 'RdBu' color scheme <br>
# with 10 levels [ref. http://colorbrewer2.org ], with the 2 centermost levels in white.

# <codecell>

# The default 'scl' map that plot.ly utilizes
def default_scl_map(x):
    """ x: some vector to be map """
    a = 1./(air.max() - air.min())     # the slope 
    return a*(x - air.min())           # the image 

N_scl = 10   # number of colors 
lwb_scl = np.min([np.floor(air.min()),-np.ceil(air.max())])  # lower bound
upb_scl = np.max([-np.floor(air.min()),np.ceil(air.max())])  # upper bound
i_scl = default_scl_map(np.linspace(lwb_scl,upb_scl,N_scl))  # scaled levels to be used
#print lwb_scl, upb_scl

# colorbrewer2.0 RdBu 
cmap = np.array([ [103, 0, 31], [178, 24, 43], [214, 96, 77], [244, 165, 130],
                  [253, 219, 199], [209, 229, 240], [146, 197, 222], 
                  [67, 147, 195], [33, 102, 172], [5, 48, 9] ])

# invert colormap so that negative temperature anomalies are in blue  
cmap = cmap[::-1,:]

# with white at 2 centermost positions
cmap[4,:] = [255,255,255]
cmap[5,:] = [255,255,255]

# now we are ready to make a plot.ly dictionary for our temperature variable
data_air = {
             'x': lon,
             'y': lat,
             'z': air, 
             'type': 'heatmap',
             'scl': # scaled level, rgb('cmap') , add ',' in-between 'cmap' elements
                    [[i_scl[i], "rgb("+','.join(map(str,cmap[i,:]))+")"] for i in range(N_scl)], 
            }
# 'text' option does not work when using the 'heatmap' type.

# <headingcell level=2>

# 5) Plot the results!

# <codecell>

# Concatenate the coastline/country list (above) and temperature list (below)
data = data_cc+[data_air]
    
py.ioff()         # turn off new tab pop-up in browser
    
# call plot.ly, save fig as 'polar_vortex'
plot_out = py.plot(
                   data,
                   layout=layout_style,
                   filename='polar_vortex',
                   fileopt='overwrite'
                  )
    
# show plot in notebook using 'plot_out'
print('\n This plot is available at '+plot_out['url'])
show_plot(plot_out['url'])

# <markdowncell>

# **Conclusion:** Cold in Eastern North America, quite warm in Alaska, Greenland and most Europe.  

