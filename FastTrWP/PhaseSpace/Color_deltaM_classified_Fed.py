''' Create plot for delta mag and color for set of transients for science justification.

Use:
(1) Type Ia blue (17cbv)
(2) Type Ia red (11fe)
(3) KN Blue (GW170817)
(4) KN Red (Model of Gw170817 w/o blue)
(5) Type IIb cooling curve (16gkg to start, add more.)
(6) FBOTs [gold templates]
(8) Type II light curve.
(9) COULD USE: just a larger sample of 'normal' events.

NOTE: for now take g and i band colors... 
Can adjust the delta t between g-band.
Should generalize to more colors etc.'''
import os
import pickle
from astropy.io import ascii
import numpy as np
from numpy.polynomial import polynomial as P
import sys
from itertools import permutations, repeat
from Calculate_ColorDelta import Calculate_ColorDelta
#Bokeh setup
from bokeh.plotting import figure,output_file,show, ColumnDataSource, gridplot
import glob
from sklearn.svm import SVC
import pandas as pd

from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn import datasets
from sklearn.svm import SVC
import pylab as pl
from sklearn.metrics import accuracy_score

import matplotlib as mpl

#from Color_deltaM_PhaseSpace import *

TRAINING = 10
#REFIT = True
REFIT = False


def plotBK(data):
    colorIa_e, deltamagIa_e,\
        colorIa2, deltamagIa2 ,\
        colorIa_l, deltamagIa_l,\
        colorIa2_l, deltamagIa2_l,\
        colorGW, deltamagGW,\
        colorGWr, deltamagGWr,\
        colorIIb_e, deltamagIIb_e,\
        colorIIb_l, deltamagIIb_l,\
        colorFB2, deltamagFB2, \
        colorFB3, deltamagFB3,\
        colorFB4, deltamagFB4,\
        colorFB3b, deltamagFB3b,\
        colorFB4b, deltamagFB4b,\
        colorFB6, deltamagFB6,\
        colorFB7, deltamagFB7,\
        colorFB5, deltamagFB5,\
        normals,nnormals\
        = data


    output_file("LSST_PhaseSpace2.html")
    tools = ['pan,box_zoom,wheel_zoom,box_select,lasso_select,crosshair,reset,save']
    
    p = figure(title="dt1 = 0min; dt2=3 hours", x_axis_label='Delta (mag)', y_axis_label='g - i',tools=tools)
    
    #Kilonovae
    p.line(deltamagGW,colorGW,legend='Kilonova (blue)',line_width=0.5,line_color='Black')
    p.circle(deltamagGW,colorGW,legend='Kilonova (blue)',line_width=0.5,size=2.5,fill_color='Blue',line_color='Black')
    p.line(deltamagGWr,colorGWr,legend='Kilonova (red)',line_width=0.5,line_color='Black')
    p.circle(deltamagGWr,colorGWr,legend='Kilonova (red)',line_width=0.5,size=2.5,fill_color='Red',line_color='Black')
    
    #Type IIb SN
    p.line(deltamagIIb_e,colorIIb_e,legend='Cooling Envelope',line_width=0.5,line_color='Black')
    p.circle(deltamagIIb_e,colorIIb_e,legend='Cooling Envelope',line_width=0.5,size=2.5,fill_color='Cyan',line_color='Black')
    
    p.line(deltamagIIb_l,colorIIb_l,legend='Cooling Envelope',line_width=0.5,line_color='Black')
    p.circle(deltamagIIb_l,colorIIb_l,legend='Cooling Envelope',line_width=0.5,size=2.5,fill_color='Cyan',line_color='Black')
    
    #Type Ia SN: Early
    p.line(deltamagIa_e,colorIa_e,legend='Type Ia (early)',line_width=0.5,line_color='Black')
    p.circle(deltamagIa_e,colorIa_e,legend='Type Ia (early)',line_width=0.5,size=2.5,fill_color='Purple',line_color='Black')
    p.line(deltamagIa2,colorIa2,legend='Type Ia (early)',line_width=0.5,line_color='Black')
    p.circle(deltamagIa2,colorIa2,legend='Type Ia (early)',line_width=0.5,size=2.5,fill_color='Purple',line_color='Black')
    
    #Type Ia SN: Normal
    p.line(deltamagIa_l,colorIa_l,legend='Type Ia (peak)',line_width=0.5,line_color='Black')
    p.circle(deltamagIa_l,colorIa_l,legend='Type Ia (peak)',line_width=0.5,size=2.5,fill_color='Grey',line_color='Black')
    p.line(deltamagIa2_l,colorIa2_l,legend='Type Ia (peak)',line_width=0.5,line_color='Black')
    p.circle(deltamagIa2_l,colorIa2_l,legend='Type Ia (peak)',line_width=0.5,size=2.5,fill_color='Grey',line_color='Black')
    
    #FBOTS:
    p.line(deltamagFB2,colorFB2,legend='FBOT rise',line_width=0.5,line_color='Black')
    p.circle(deltamagFB2,colorFB2,legend='FBOT rise',line_width=0.5,size=2.5,fill_color='Yellow',line_color='Black')
    p.line(deltamagFB3,colorFB3,legend='FBOT rise',line_width=0.5,line_color='Black')
    p.circle(deltamagFB3,colorFB3,legend='FBOT rise',line_width=0.5,size=2.5,fill_color='Yellow',line_color='Black')
    p.line(deltamagFB4,colorFB4,legend='FBOT rise',line_width=0.5,line_color='Black')
    p.circle(deltamagFB4,colorFB4,legend='FBOT rise',line_width=0.5,size=2.5,fill_color='Yellow',line_color='Black')
    
    
    p.line(deltamagFB3b,colorFB3b,legend='FBOT decline',line_width=0.5,line_color='Black')
    p.circle(deltamagFB3b,colorFB3b,legend='FBOT decline',line_width=0.5,size=2.5,fill_color='IndianRed',line_color='Black')
    p.line(deltamagFB4b,colorFB4b,legend='FBOT decline',line_width=0.5,line_color='Black')
    p.circle(deltamagFB4b,colorFB4b,legend='FBOT decline',line_width=0.5,size=2.5,fill_color='IndianRed',line_color='Black')
    
    p.line(deltamagFB6,colorFB6,legend='FBOT decline',line_width=0.5,line_color='Black')
    p.circle(deltamagFB6,colorFB6,legend='FBOT decline',line_width=0.5,size=2.5,fill_color='IndianRed',line_color='Black')
    p.line(deltamagFB7,colorFB7,legend='FBOT decline',line_width=0.5,line_color='Black')
    p.circle(deltamagFB7,colorFB7,legend='FBOT decline',line_width=0.5,size=2.5,fill_color='IndianRed',line_color='Black')
    
    #Rapidly declining Type I SN
    p.line(deltamagFB5,colorFB5,legend='Fast decliner',line_width=0.5,line_color='Black')
    p.circle(deltamagFB5,colorFB5,legend='Fast decliner',line_width=0.5,size=2.5,fill_color='Salmon',line_color='Black')
    for sn in normals:
        p.circle(sn[0], sn[1], line_width=0.5,size=2.5,fill_color=None,line_color='Black')
    
    
    show(p)

def readdata():

    rawdata = {}
    rawdata['annotations'] = {}
    
    #-KILONOVAE------------------------------------------------------------
    rawdata['GW'] = ascii.read('LightCurves/Formatted/GW170817.dat')
    rawdata['annotations']['GW'] = 'kilonova'
    #colorGW,deltamagGW = Calculate_ColorDelta(dataGW['g'],dataGW['i'],delta_t,color_t)

    #-RED-KILONOVAE-----------------------------------------------
    rawdata['GWr'] = ascii.read('LightCurves/Formatted/GW170817red.dat')
    rawdata['annotations']['GWr'] = 'red kilonova'
    #colorGWr,deltamagGWr = Calculate_ColorDelta(dataGWr['g'],dataGWr['i'],delta_t,color_t)


    #-RISING COOLING ENVELOPE (16gkg)------------------------------
    rawdata['IIb'] = ascii.read('LightCurves/Formatted/SN2016gkg_poly_e.dat')
    rawdata['annotations']['IIb'] = 'RISING COOLING ENVELOPE 16gkg'
    #colorIIb_e,deltamagIIb_e = Calculate_ColorDelta(dataIIb['g'],dataIIb['i'],delta_t,color_t)

    #-DECLINING COOLING ENVELOPE (16gkg)------------------------------
    rawdata['IIb_l'] = ascii.read('LightCurves/Formatted/SN2016gkg_poly_l.dat')
    rawdata['annotations']['IIb_l'] = 'DECLINING COOLING ENVELOPE 16gkg'
    #colorIIb_l,deltamagIIb_l = Calculate_ColorDelta(dataIIb_l['g'],dataIIb_l['i'],delta_t,color_t)


    #-EARLY TYPE IA SN--(< 6 days from explosion)------------------
    #-SN2017cbv----------------------------------------------------
    rawdata['Ia'] = ascii.read('LightCurves/Formatted/SN2017cbv_poly_e.dat')
    rawdata['annotations']['Ia'] = 'EARLY TYPE IA SN--(< 6 days from explosion)'
    #colorIa_e,deltamagIa_e = Calculate_ColorDelta(dataIa['g'],dataIa['i'],delta_t,color_t)

    #-SN2011fe------------------------------------------------------
    rawdata['Ia2'] = ascii.read('LightCurves/Formatted/SN2011fe_poly_e.dat')
    rawdata['annotations']['Ia2'] = 'sn2011fe'
    #colorIa2,deltamagIa2 = Calculate_ColorDelta(dataIa2['g'],dataIa2['i'],delta_t,color_t)


    #-NORMAL TYPE IA SN (> 10 days post-explosion-----------------_
    #-SN2017cbv----------------------------------------------------
    rawdata['Ia_l'] = ascii.read('LightCurves/Formatted/SN2017cbv_poly_l.dat')
    rawdata['annotations']['Ia_l'] = 'NORMAL TYPE IA SN (> 10 days post-explosion 17cbv'
    #colorIa_l,deltamagIa_l = Calculate_ColorDelta(dataIa_l['g'],dataIa_l['i'],delta_t,color_t)

    #-SN2011fe------------------------------------------------------
    rawdata['Ia2_l'] = ascii.read('LightCurves/Formatted/SN2011fe_poly_l.dat')
    rawdata['annotations']['Ia2_l'] = '2011fe ?? what is the difference between this and Ia2??'
    #colorIa2_l,deltamagIa2_l = Calculate_ColorDelta(dataIa2_l['g'],dataIa2_l['i'],delta_t,color_t)


    #-FBOTS: RISING PHASE-----------------------------------------
    #-PS10bjp-----------------------------------------------------
    rawdata['FB2'] = ascii.read('LightCurves/Formatted/PS10bjp_rise.dat')
    rawdata['annotations']['FB2'] = 'PS10bjp fbot'
    #colorFB2,deltamagFB2 = Calculate_ColorDelta(dataFB2['g'],dataFB2['i'],delta_t,color_t)

    #-PS10ah-----------------------------------------------------
    rawdata['FB3'] = ascii.read('LightCurves/Formatted/PS10ah_rise.dat')
    rawdata['annotations']['FB3'] = 'PS10ah fbot'
    #colorFB3,deltamagFB3 = Calculate_ColorDelta(dataFB3['g'],dataFB3['i'],delta_t,color_t)

    #-PS12brf-----------------------------------------------------
    rawdata['FB4'] = ascii.read('LightCurves/Formatted/PS12brf_rise.dat')
    rawdata['annotations']['FB4'] = ' PS12br fbot'
    #colorFB4,deltamagFB4 = Calculate_ColorDelta(dataFB4['g'],dataFB4['i'],delta_t,color_t)


    #-FBOTS: DECLINING PHASE--------------------------------------
    #-PS10ah------------------------------------------------------
    rawdata['FB3b'] = ascii.read('LightCurves/Formatted/PS10ah_decline.dat')
    rawdata['annotations']['FB3b'] = 'PS10ah declining'
    #colorFB3b,deltamagFB3b = Calculate_ColorDelta(dataFB3b['g'],dataFB3b['i'],delta_t,color_t)

    #-PS12brf------------------------------------------------------
    rawdata['FB4b'] = ascii.read('LightCurves/Formatted/PS12brf_decline.dat')
    rawdata['annotations']['FB4b'] = 'PS12brf fbot'
    #colorFB4b,deltamagFB4b = Calculate_ColorDelta(dataFB4b['g'],dataFB4b['i'],delta_t,color_t)

    #-PS11qr--------------------------------------------------------
    rawdata['FB6'] = ascii.read('LightCurves/Formatted/PS11qr_decline.dat')
    rawdata['annotations']['FB6'] = 'PS11qr fbot'
    #colorFB6,deltamagFB6 = Calculate_ColorDelta(dataFB6['g'],dataFB6['i'],delta_t,color_t)

    #-AT2018cow------------------------------------------------------
    rawdata['FB7'] = ascii.read('LightCurves/Formatted/AT2018cow_decline.dat')
    rawdata['annotations']['FB7'] = 'AT2018cow fbot'
    #colorFB7,deltamagFB7 = Calculate_ColorDelta(dataFB7['g'],dataFB7['i'],delta_t,color_t)

    
    rawdata['FB5'] = ascii.read('LightCurves/Formatted/SN2005ek_decline.dat')
    rawdata['annotations']['FB5'] = 'SN2005ek decline'


    #-------------------------NORMALS--------------------------------
    #----------------------------------------------------------------

    #-IIs------------------------------------------------------------


    
    for i,f in enumerate(glob.glob("IIs/*.dat")):
        rawdata['sn_II_%i'%i] = ascii.read(f, format='csv', delimiter="\t")  
        rawdata['annotations']['snII_%i'%i] = f

    #-Ia's------------------------------------------------------------    
    for i,f in enumerate(glob.glob("IAs/*.csv")):
        rawdata['sn_%i'%i] = ascii.read(f)
        rawdata['annotations']['sn_%i'%i] = f
        
    #returnvalue = (dataGW, dataGWr, dataIIb, dataIIb_l,
     #              dataIa, dataIa2, dataIa_l, dataIa2_l,
     #              dataFB2, dataFB3, dataFB4, dataFB3b,
     #              dataFB4b, dataFB6, dataFB7, dataFB5,
     #              [sn for sn in sne])
    return rawdata


    
def getcdt(datain,
           delta_t=3.0, color_t=.5):
    
    '''dataGW, dataGWr,dataIIb, dataIIb_l, \
        dataIa, dataIa2, dataIa_l, dataIa2_l,\
        dataFB2, dataFB3, dataFB4, dataFB3b, \
        dataFB4b, dataFB6, dataFB7, dataFB5, \
        sne = datain
    '''
           #Set some initial parameters:
    #delta_t = 3.0 #set in terms of number of hours. Should always be in units of 0.5 hours.
    #color_t = 0. #set in hours. Must either be 0 or 0.5 (for 30 minutes)

    ######Reading in all the files: This will swap out when I finalize all of the new template models.
    ######I am almost done writing out all of the new templates so we won't do the interpolation in code,
    ######Instead it will just read in the templates whiles have been interpolated to 0.5 days.

    ######This version of the code contains: blue and red KN. 16gkg cooling envelope, Type Ia with blue bump
    ######The Type Ia is separated into (a) early LC (first 5 days post explosion and (b) Standard (10+ days post).

    
    returnvalue = []

    Calculate_ColorDelta(datain['GW']['g'], datain['GW']['i'],
                         delta_t, color_t)
    for k in datain.keys():
        returnvalue.append(Calculate_ColorDelta(datain[k]['g'],
                                                datain[k]['i'],
                                                delta_t,
                                                color_t))
    return returnvalue






dt1 = [1.5]#, 0.5, 3.5, 4.5, 6.5]#[0.5,1,1.5,2.5,3.5,4.5,5.5,6.5] #gap between obs in same filter
dt2 = [0.5]#, 0, 1, 2] #gap between filters
data = readdata()

#print(data['GWr'])
newdata = {}

for i,did in enumerate(data):
    if did == 'annotations':
        continue
    #[:-2]):

    d = data[did]   
    indices = d['time-rel'] < 15*24
    indices = d['time-rel'] > -21*24    
    newdata[did] = d[indices]

    #pl.plot(d[d['time-rel'])
#newdata.append(data[-2])
#newdata.append(data[-1])               
               
    #d = d[d['time-rel'] < 30]
#pl.ion()

if __name__ == '__main__':
    kernel = 1.0 * RBF([1.0, 1.0])
    if REFIT:
        clf = GaussianProcessClassifier(kernel)
    #print(np.array([k for k in newdata.keys()]))
    nspecial = np.sum(~np.array([k.startswith('sn_') for k in newdata.keys()]))
    #number of normal IIs
    nII = np.sum(np.array([k.startswith('sn_II_') for k in newdata.keys()]))
    for i in [(t1, t2) for t1 in dt1 for t2 in dt2]:
        print(i)
        thisdir = "GPclassifier_%.1f_%.1f"%i
        os.system("mkdir " + thisdir)
        tmp2 = getcdt(newdata,
                  delta_t=i[0], color_t=i[1])
        
        
        color = np.hstack([t[0] for t in tmp2[:nspecial]])
        shape = np.hstack([t[1] for t in tmp2[:nspecial]])
        #nsnIa = np.array([len(tmp2[i-2]) for i in range(4)]).sum()
        iicolor = np.hstack([t[0] for t in tmp2[nspecial:nspecial+nII]])
        iishape = np.hstack([t[1] for t in tmp2[nspecial:nspecial+nII]])

        N = color.shape[0]
        nii = iicolor.shape[0]

        iacolor = np.hstack([t[0] for t in tmp2[nspecial+nII:]])
        iashape = np.hstack([t[1] for t in tmp2[nspecial+nII:]])
        nia = iacolor.shape[0]        

        np.random.seed(123)
        
        # pick boring transients: 1/3 II and 1/2 Ia
        randIa = np.random.randint(0, len(iacolor), int(2 * N / 3))
        randII = np.random.randint(0, len(iicolor), int(N / 3))

        plotcolor = ['IndianRed'] * N + ['c'] * N
        label = np.hstack([np.zeros(N), np.ones(N)])
        #print(N, len(label))
        
        color = np.hstack([color, iicolor[randII], iacolor[randIa]])
        shape = np.hstack([shape, iishape[randII], iashape[randIa]])
        #print(color.min(), color.max())
        #print(shape.min(), shape.max())        
                # pick a TRAINING fraction 
        indx = np.random.randint(0, len(label), size=int(N/TRAINING))
        #print(min(indx), max(indx))

        phasespace_complete = np.array([shape, color]).T
        #print(phasespace_complete.shape)
        
        X = phasespace_complete[indx] #training set
        y = label[indx] #training labels


        fig = pl.figure(figsize=(10,5))
        ax = fig.add_subplot(111)

        ax.plot(shape[-N:-N + int(N / 3)], color[-N:-N + int(N / 3)],
                '.', c='mediumpurple',
                alpha=0.05, ms=2)
        ax.plot(shape[-N + int(N / 3):], color[-N + int(N / 3):], '.', c='SteelBlue',
                alpha=0.03, ms=2)            
        ax.plot(shape[:-N], color[:-N], '.', c='IndianRed', alpha=0.5, ms=2)
        ax.set_xlim(-0.7, 1.45)
        ax.set_ylim(-1.2, 4.2)
        ax.text(1.0, 3.0, r"$\Delta t_1 = %.1f$"%i[0], fontsize=18)
        ax.text(1.0, 2.5, r"$\Delta t_2 = %.1f$"%i[1], fontsize=18)
        #pl.draw()
        #nsnIa = np.array([len(tmp2[i-2]) for i in range(4)]).sum()
        
        #label = np.zeros(len(color))
        #label[:nsnIa] = 1

        if REFIT:
            clf.fit(X, y)
            pickle.dump(clf, open(thisdir+"/GPmodel.pkl", 'wb'))
        else:
            # load the model from disk
            clf = pickle.load(open(thisdir+"/GPmodel.pkl", 'rb'))
            

        y_pred = clf.predict(X)
        #in sample accuracy
        accuracy = accuracy_score(y, y_pred)
        print("Accuracy (train) for %0.1f%% " % (accuracy * 100))

        #calculate probability everywhere in the phase space
        '''
        xx = np.linspace(phasespace_complete[:,0].min()-0.1,
                         phasespace_complete[:,0].max()+0.1, 100)
        yy = np.linspace(phasespace_complete[:,1].min()-0.1,
                         phasespace_complete[:,1].max()+0.1, 100).T
        '''
        resolution = 50
        xx = np.linspace(-0.7, 1.45, resolution)
        yy = np.linspace(-1.2, 4.2, resolution)
        xx, yy = np.meshgrid(xx, yy)
        Xfull = np.c_[xx.ravel(), yy.ravel()]
        probas = clf.predict_proba(Xfull)

        # View probabilities:
        n_classes = np.unique(y_pred).size
        
        imshow_handle = ax.imshow(probas[:, 0].reshape((resolution, resolution)),
                                   extent=(-0.7, 1.45,
                                           -1.2, 4.2),
                                  norm=mpl.colors.Normalize(vmin=0, vmax=1),
                                  cmap='BrBG', aspect='auto', origin='lower')

        
        idx = (y_pred == 0)
        if idx.any():
            ax.scatter(X[:, 0], X[:, 1], marker='o', c='none', edgecolor='k',
                       alpha=0.5)

        ax.set_xlabel(r"$\Delta$ mag", fontsize=18)
        ax.set_ylabel(r"color", fontsize=18)
        cbaxes = fig.add_axes([0.125, 0.9, 0.775, 0.02]) 
        cbaxes.set_title("Probability of unusual transient")


        cb = pl.colorbar(imshow_handle, orientation='horizontal', cax=cbaxes,
                         ticks=[0, 0.25, 0.5, 0.75, 1],
                         norm=mpl.colors.Normalize(vmin=0.01, vmax=0.99))
        cb.ax.xaxis.set_ticks_position('top')
        cb.ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '1'])
        
        #cb.ax.set_xlim(0,1)
        #, boundaries=[0,1])
        #pl.draw()
        mask = np.array([~(i in indx)
                            for i in range(len(np.array([shape, color]).T))])

        y_predOUT = clf.predict(np.array([shape, color]).T[mask])
        yOUT = label[mask]
        print (len(yOUT),  accuracy_score(yOUT, y_predOUT))
        pl.savefig(thisdir + "/GPclassifier_%.1f_%.1f.png"%i)
        
        '''
        XX, YY = np.mgrid[shape.min():shape.max():200j,
                      color.min():color.max():200j]

        clf = SVC(gamma='auto', kernel='rbf')
        clf.fit(np.array(list(zip(np.log10(shape + 1), color))), label)
        Z = clf.decision_function(np.c_[np.log10(XX.ravel() + 1), YY.ravel()])
        Z = Z.reshape(XX.shape)
        pl.contour(XX, YY, Z, colors=['k', 'k', 'k'],
        linestyles=['--', '-', '--'], levels=[-.5, 0, .5])

        #pl.pcolormesh(XX, YY, Z > 0, cmap=pl.cm.bone)
        ## Circle out the test dat
        #pl.plot(shape[nsnIa:], color[nsnIa:],'.', c="SteelBlue")
        pl.title(r"dt_1=%.1f dt_2=%.1f"%(i[0],i[1]))
        pl.savefig("dt_1=%.1fdt_2=%.1f.png"%(i[0],i[1]))
    pl.show()

        
   for k in range(n_classes):
        plt.title("Class %d" % k)
        
        imshow_handle = plt.imshow(probas[:, k].reshape((100, 100)),
                                   extent=(X.T[0].min(),
            X.T[0].max(), X.T[1].min(), X.T[1].max()), 
            aspect='auto', origin='lower')
        plt.xticks(())
        plt.yticks(())
        idx = (y_pred == k)
        if idx.any():
            plt.scatter(X[idx, 0], X[idx, 1], marker='o', c='w', edgecolor='k')

        
        imshow_handle = plt.imshow(probas[:, k].reshape((100, 100)),
                                   extent=(X.T[0].min(),
            X.T[0].max(), X.T[1].min(), X.T[1].max()), 
            aspect='auto', origin='lower')
        plt.xticks(())
        plt.yticks(())
        idx = (y_pred == k)
        if idx.any():
            plt.scatter(X[idx, 0], X[idx, 1], marker='o', c='w', edgecolor='k')

plt.title("Probability")
Xfull = np.c_[xx.ravel(), yy.ravel()]
if 1:
    y_pred = clf.predict(X)
    accuracy = accuracy_score(y, y_pred)
    print("Accuracy (trai) for  %0.1f%% " %( accuracy * 100))

    # View probabilities:
    probas = clf.predict_proba(Xfull)
    n_classes = np.unique(y_pred).size
    for k in range(n_classes):
        plt.title("Class %d" % k)
        
        imshow_handle = plt.imshow(probas[:, k].reshape((100, 100)),
                                   extent=(X.T[0].min(),
            X.T[0].max(), X.T[1].min(), X.T[1].max()), 
            aspect='auto', origin='lower')
        plt.xticks(())
        plt.yticks(())
        idx = (y_pred == k)
        if idx.any():
            plt.scatter(X[idx, 0], X[idx, 1], marker='o', c='w', edgecolor='k')

plt.title("Probability")
if 1:
'''
