''' Makes Figure 4 and 5 or PrestoColor

Complete code including model fitting in Color_deltaM_clalssify.py
'''

import os
import pickle
import astropy
from astropy import table
from astropy.io import ascii
import numpy as np
from numpy.polynomial import polynomial as P
import sys
from itertools import permutations, repeat
from Calculate_ColorDelta_Fed import Calculate_ColorDelta
import glob
import pandas as pd

from sklearn.preprocessing import Imputer 
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn import datasets
from sklearn.svm import SVC
import pylab as pl
from sklearn.metrics import accuracy_score
from Calculate_ColorDelta_Fed import Calculate_ColorDelta
from utils import readdata
from utils import getcdt
from utils import dropna
from utils import TMIN, TMAX, SEED, TRAINING, uncertainties, prepdata

from bokeh.layouts import column
from bokeh.models import CustomJS, ColumnDataSource, Slider, Select
from bokeh.plotting import figure, output_file, show
from bokeh.models.glyphs import VBar
import pandas as pd
from bokeh.layouts import row, column
from bokeh.models import CustomJS, Slider
from bokeh.plotting import Figure
from bokeh.models.widgets.buttons import Button
from bokeh.themes import built_in_themes
from bokeh.io import curdoc
from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead, Label

from bokeh.plotting import figure, output_file, show, ColumnDataSource
import matplotlib as mpl
pl.style.use('dark_background')

mpl.rcParams.update({'font.size': 18})
#pl.ion()
#from Color_deltaM_PhaseSpace import *
PLASTICC = False
filt2 = {'g':'i', 'r':'i'}

xmin, xmax = -0.6, 1.05
ymin, ymax = -1.05, 3.5
stepc = 0.5
steps = 1

r1 = np.arange(0.0, 4.5, stepc)#0.5)
r2 = np.arange(1., 9., steps)#1)#np.ones_like(np.arange(1, 8, 1))

bestvalues = (0.5, 4.0)
N = 0        

assert(len(sys.argv) == 2), '''


proper usage:
python paperInteractFig.py g 
or 
python paperInteractFig.py r '''

ofilt = sys.argv[1]

def getuncertainties(filt):
    unc = []
    for t in r1:
        unc.append([ uncertainties[filt][t][20],
                     uncertainties[filt][t][22],
                     uncertainties[filt][t][24]])
        
    return unc

def classifierfig(filt='g'):
    global N
    data = readdata(PLASTICC)
    newdata = prepdata(data, TMIN, TMAX)
    

    plots = [] #will host all data
    ids = [] #id to change data source
    nspecial = np.sum(~np.array([k.startswith('sn_')
                                 for k in newdata.keys()])) 
    print ("special lightcurves", nspecial)
    #number of normal IIs
    nII = np.sum(np.array([k.startswith('sn_II_') for k in newdata.keys()]))

    for j,tmp0 in enumerate(r1):
     for k,tmp1 in enumerate(r2):
        i = k + j * len(r2)
        dts = (tmp0, tmp1)
        
        thisdir = "noplasticc/GPclassifier_%si"%filt + \
                  "/GPclassifier_%.1f_%.1f"%(dts[1], dts[0])
        tmp2 = getcdt(newdata, delta_t=dts[1],
                      color_t=dts[0], filt=filt)
        
        color = np.hstack([t[0] for t in tmp2[:nspecial]])
        shape = np.hstack([t[1] for t in tmp2[:nspecial]])
        #cleannan
        mask = np.isnan(color) * np.isnan(shape)
        color = color[~mask]
        shape = shape[~mask]

        #nsnIa = np.array([len(tmp2[i-2]) for i in range(4)]).sum()
        iicolor = np.hstack([t[0] for t in tmp2[nspecial:nspecial+nII]])
        iishape = np.hstack([t[1] for t in tmp2[nspecial:nspecial+nII]])
        #cleannan
        mask = np.isnan(iicolor) * np.isnan(iishape)
        iicolor = iicolor[~mask]
        iishape = iishape[~mask]        
        nii = iicolor.shape[0]

        iacolor = np.hstack([t[0] for t in tmp2[nspecial+nII:]])
        iashape = np.hstack([t[1] for t in tmp2[nspecial+nII:]])
        #cleannan
        mask = (np.isnan(iacolor) |  np.isnan(iashape))
        iacolor = iacolor[~mask]
        iashape = iashape[~mask]        
        nia = iacolor.shape[0]        
        np.random.seed(SEED)
        N = color.shape[0]
        
        # pick boring transients: 1/3 II and 1/2 Ia
        randIa = np.random.randint(0, len(iacolor), int(2 * N / 3))
        randII = np.random.randint(0, len(iicolor), int(N / 3) + 1)

        plotcolor = ['IndianRed'] * N + ['c'] * N
        label = np.hstack([np.zeros(N), np.ones(N)])
        #print(N, len(label))
        
        color = list(np.hstack([color, iicolor[randII], iacolor[randIa]]))
        shape = list(np.hstack([shape, iishape[randII], iashape[randIa]]))
        # pick a TRAINING fraction
        plots.append({'x':shape, 'y':color,
                      'idt':"%.1f-%.1f"%(dts[0], dts[1])})
        
        #plots['shape'] = plots['shape'] + shape
        #plots['c'] = plots['c'] + ([dts[0]] * len(color))
        #plots['s'] = plots['s'] + ([dts[1]] * len(color))  
        ids.append("%.1f-%.1f"%(dts[0], dts[1]))
        #alldata.append()
        phasespace_complete = np.array([shape, color]).T
    return plots, ids, N

plots, ids, N = classifierfig(filt=ofilt)

u = getuncertainties(ofilt)

#plots['shape'] = np.array(plots['shape'])
#plots['color'] = np.array(plots['color'])
#plots['s'] = np.array(plots['s'])
#plots['c'] = np.array(plots['c'])

#  if p['idt'] in "0.5_4.0" ])
i = [i for i,p in enumerate(plots) if p['idt'] in "0.5-4.0"][0]
j = np.where(r2 == 4.)[0][0]

colors = ["#%02x%02x%02x" % (205, 120, 50)] * N + \
    ["#%02x%02x%02x" % (170, 170, 195)] * N 

source = ColumnDataSource(data=dict(x=plots[i]['x'],
                                    y=plots[i]['y'],
                                    colors=colors,
                                    alpha=[0.3] * N + [0.1] * N))
u20 = ColumnDataSource(data=dict(
    x=[-u[j][0], u[j][0]],
    y=[-1,-1]))
u22 = ColumnDataSource(data=dict(
    x=[-u[j][1], u[j][1]],
    y=[-1.25,-1.25]))
u24 = ColumnDataSource(data=dict(
    x=[-u[j][2], u[j][2]],
    y=[-1.5,-1.5]))

#alldata = ColumnDataSource(data=plots)
#dict(x=plots['shape'],
#                                     c = plots['c'],
#                                     s = plots['s'],
#                                     y=plots['color']#[(plots['c'] == 0.5) *
#                                                    #(plots['s'] == 4)]

curdoc().theme = 'dark_minimal'

plot = figure(y_range=(-1.7, 3.5), x_range=(-0.55, 1.05),
              plot_width=400, plot_height=400)


plot.xaxis.axis_label = 'Δ %s'%ofilt
plot.yaxis.axis_label = '%s-i'%ofilt

plot.scatter('x', 'y', fill_color="colors", fill_alpha="alpha",
             line_color=None,
             source=source)

plot.line(x='x', y='y', line_color="White", line_cap='butt',
             source=u20)
plot.line(x=[-u20.data['x'][0], -u20.data['x'][0]],
          y=[u20.data['y'][0]-0.05, u20.data['y'][0]+0.05], line_color="White")
plot.line(x=[u20.data['x'][0], u20.data['x'][0]],
          y=[u20.data['y'][0]-0.05, u20.data['y'][0]+0.05], line_color="White")
plot.add_layout(Label(x=-0.5, y=-1.05, text='ε(%s=20)'%ofilt, text_color="White",
                      text_font_size='8pt'))

plot.line(x='x', y='y', line_color="White", line_cap='butt',
             source=u22)
plot.line(x=[-u22.data['x'][0], -u22.data['x'][0]],
          y=[u22.data['y'][0]-0.05, u22.data['y'][0]+0.05], line_color="White")
plot.line(x=[u22.data['x'][0], u22.data['x'][0]],
          y=[u22.data['y'][0]-0.05, u22.data['y'][0]+0.05], line_color="White")
plot.add_layout(Label(x=-0.5, y=-1.3, text='ε(%s=22)'%ofilt, text_color="White",
                      text_font_size='8pt'))

plot.line(x='x', y='y', line_color="White", line_cap='butt',
             source=u24)
plot.line(x=[-u24.data['x'][0], -u24.data['x'][0]],
          y=[u24.data['y'][0]-0.05, u24.data['y'][0]+0.05], line_color="White")
plot.line(x=[u24.data['x'][0], u24.data['x'][0]],
          y=[u24.data['y'][0]-0.05, u24.data['y'][0]+0.05], line_color="White")
plot.add_layout(Label(x=-0.5, y=-1.55, text='ε(%s=24)'%ofilt, text_color="White",
                      text_font_size='8pt'))

plot.add_layout(Arrow(end=OpenHead(line_color="White", size=10),
                      x_start=-0.350, y_start=3.3, x_end=-0.45, y_end=3.3,
                      line_color="White"))

plot.add_layout(Arrow(end=OpenHead(line_color="White", size=10),
                      x_start=0.8, y_start=3.3, x_end=0.9, y_end=3.3,
                      line_color="White"))

plot.add_layout(Arrow(end=OpenHead(line_color="firebrick", size=10),
                      x_start=0.9, y_start=2.5, x_end=0.9, y_end=3.0,
                      line_color="firebrick"))

plot.add_layout(Arrow(end=OpenHead(line_color="steelblue", size=10),
                      x_start=0.9, y_start=-0.7, x_end=0.9, y_end=-1.4,
                      line_color="steelblue"))


plot.add_layout(Label(x=-0.3, y=3.2, text='decline', text_color="White"))
plot.add_layout(Label(x=0.65, y=3.2, text='rise', text_color="White"))

plot.add_layout(Label(x=0.93, y=2.05, text='red', angle=90, angle_units="deg",
                      text_color="firebrick"))
plot.add_layout(Label(x=0.93, y=-0.55, text='blue', angle=90, angle_units="deg",
                      text_color="steelblue"))


####WIDGETS

cslider = Slider(start=r1[0], end=4, value=bestvalues[0],
                 step=stepc, title="Tcolor")
sslider = Slider(start=r2[0], end=8, value=bestvalues[1],
                 step=steps, title="Tshape")

phase_slider = Slider(start=0, end=6.4, value=0, step=.1, title="Phase")
offset_slider = Slider(start=-5, end=5, value=0, step=.1, title="Offset")

callback = CustomJS(args=dict(source=source, c=cslider, s=sslider,
                              cc=plots, ids=ids),
                    code='''
    var data = source.data;
    var allplots = cc;
    var ids=ids;
    var csl = c.value;
    var ssl = s.value;
    console.log(ids);
    var tmp = csl.toFixed(1) + "-" + ssl.toFixed(1);

    console.log(tmp, ids.indexOf(tmp), allplots[ids.indexOf(tmp)]);
    data['x'] =  allplots[ids.indexOf(tmp)]['x'];
    data['y'] =  allplots[ids.indexOf(tmp)]['y'];
    source.change.emit();
''')

#cslider.js_on_change('value', callback)
sslider.js_on_change('value', callback)
cslider.js_on_change('value', callback)

# Set up layouts and add to document

output_file("%si.html"%ofilt, title="slider.py example")

layout = column(
    plot, cslider, sslider)

show(layout)


