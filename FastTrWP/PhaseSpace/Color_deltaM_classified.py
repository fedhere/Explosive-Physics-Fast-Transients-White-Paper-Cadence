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
from Color_deltaM_PhaseSpace import *

def readdata():
    #-KILONOVAE------------------------------------------------------------
    dataGW = ascii.read('LightCurves/Formatted/GW170817.dat')
    #colorGW,deltamagGW = Calculate_ColorDelta(dataGW['g'],dataGW['i'],delta_t,color_t)

    #-RED-KILONOVAE-----------------------------------------------
    dataGWr = ascii.read('LightCurves/Formatted/GW170817red.dat')
    #colorGWr,deltamagGWr = Calculate_ColorDelta(dataGWr['g'],dataGWr['i'],delta_t,color_t)


    #-RISING COOLING ENVELOPE (16gkg)------------------------------
    dataIIb = ascii.read('LightCurves/Formatted/SN2016gkg_poly_e.dat')
    #colorIIb_e,deltamagIIb_e = Calculate_ColorDelta(dataIIb['g'],dataIIb['i'],delta_t,color_t)

    #-DECLINING COOLING ENVELOPE (16gkg)------------------------------
    dataIIb_l = ascii.read('LightCurves/Formatted/SN2016gkg_poly_l.dat')
    #colorIIb_l,deltamagIIb_l = Calculate_ColorDelta(dataIIb_l['g'],dataIIb_l['i'],delta_t,color_t)


    #-EARLY TYPE IA SN--(< 6 days from explosion)------------------
    #-SN2017cbv----------------------------------------------------
    dataIa = ascii.read('LightCurves/Formatted/SN2017cbv_poly_e.dat')
    #colorIa_e,deltamagIa_e = Calculate_ColorDelta(dataIa['g'],dataIa['i'],delta_t,color_t)

    #-SN2011fe------------------------------------------------------
    dataIa2 = ascii.read('LightCurves/Formatted/SN2011fe_poly_e.dat')
    #colorIa2,deltamagIa2 = Calculate_ColorDelta(dataIa2['g'],dataIa2['i'],delta_t,color_t)


    #-NORMAL TYPE IA SN (> 10 days post-explosion-----------------_
    #-SN2017cbv----------------------------------------------------
    dataIa_l = ascii.read('LightCurves/Formatted/SN2017cbv_poly_l.dat')
    #colorIa_l,deltamagIa_l = Calculate_ColorDelta(dataIa_l['g'],dataIa_l['i'],delta_t,color_t)

    #-SN2011fe------------------------------------------------------
    dataIa2_l = ascii.read('LightCurves/Formatted/SN2011fe_poly_l.dat')
    #colorIa2_l,deltamagIa2_l = Calculate_ColorDelta(dataIa2_l['g'],dataIa2_l['i'],delta_t,color_t)


    #-FBOTS: RISING PHASE-----------------------------------------
    #-PS10bjp-----------------------------------------------------
    dataFB2 = ascii.read('LightCurves/Formatted/PS10bjp_rise.dat')
    #colorFB2,deltamagFB2 = Calculate_ColorDelta(dataFB2['g'],dataFB2['i'],delta_t,color_t)

    #-PS10ah-----------------------------------------------------
    dataFB3 = ascii.read('LightCurves/Formatted/PS10ah_rise.dat')
    #colorFB3,deltamagFB3 = Calculate_ColorDelta(dataFB3['g'],dataFB3['i'],delta_t,color_t)

    #-PS12brf-----------------------------------------------------
    dataFB4 = ascii.read('LightCurves/Formatted/PS12brf_rise.dat')
    #colorFB4,deltamagFB4 = Calculate_ColorDelta(dataFB4['g'],dataFB4['i'],delta_t,color_t)


    #-FBOTS: DECLINING PHASE--------------------------------------
    #-PS10ah------------------------------------------------------
    dataFB3b = ascii.read('LightCurves/Formatted/PS10ah_decline.dat')
    #colorFB3b,deltamagFB3b = Calculate_ColorDelta(dataFB3b['g'],dataFB3b['i'],delta_t,color_t)

    #-PS12brf------------------------------------------------------
    dataFB4b = ascii.read('LightCurves/Formatted/PS12brf_decline.dat')
    #colorFB4b,deltamagFB4b = Calculate_ColorDelta(dataFB4b['g'],dataFB4b['i'],delta_t,color_t)

    #-PS11qr--------------------------------------------------------
    dataFB6 = ascii.read('LightCurves/Formatted/PS11qr_decline.dat')
    #colorFB6,deltamagFB6 = Calculate_ColorDelta(dataFB6['g'],dataFB6['i'],delta_t,color_t)

    #-AT2018cow------------------------------------------------------
    dataFB7 = ascii.read('LightCurves/Formatted/AT2018cow_decline.dat')
    #colorFB7,deltamagFB7 = Calculate_ColorDelta(dataFB7['g'],dataFB7['i'],delta_t,color_t)
    dataFB5 = ascii.read('LightCurves/Formatted/SN2005ek_decline.dat')

    sne = []
    for i,f in enumerate(glob.glob("IAs/*.csv")):
        sne.append(ascii.read(f))

    returnvalue = (dataGW, dataGWr, dataIIb, dataIIb_l,
                   dataIa, dataIa2, dataIa_l, dataIa2_l,
                   dataFB2, dataFB3, dataFB4, dataFB3b,
                   dataFB4b, dataFB6, dataFB7, dataFB5,
                   [sn for sn in sne])
    return returnvalue



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
    
def getcdt(datain,
           delta_t=3.0, color_t=.5):
    
    dataGW, dataGWr,dataIIb, dataIIb_l, \
        dataIa, dataIa2, dataIa_l, dataIa2_l,\
        dataFB2, dataFB3, dataFB4, dataFB3b, \
        dataFB4b, dataFB6, dataFB7, dataFB5, \
        sne = datain
    
           #Set some initial parameters:
    #delta_t = 3.0 #set in terms of number of hours. Should always be in units of 0.5 hours.
    #color_t = 0. #set in hours. Must either be 0 or 0.5 (for 30 minutes)

    step = np.int(delta_t/0.5)
    color_step = np.int(color_t/0.5)


    ######Reading in all the files: This will swap out when I finalize all of the new template models.
    ######I am almost done writing out all of the new templates so we won't do the interpolation in code,
    ######Instead it will just read in the templates whiles have been interpolated to 0.5 days.

    ######This version of the code contains: blue and red KN. 16gkg cooling envelope, Type Ia with blue bump
    ######The Type Ia is separated into (a) early LC (first 5 days post explosion and (b) Standard (10+ days post).

    #------------------------------------------------

    color_t = 0. #set in hours. Must either be 0 or 0.5 (for 30 minutes)
    
    #----------------------------------------------------------------------
    #-READ IN TEMPLATES----------------------------------------------------
    #----------------------------------------------------------------------
    
    #-KILONOVAE------------------------------------------------------------
    colorGW,deltamagGW = Calculate_ColorDelta(dataGW['g'],dataGW['i'],delta_t,color_t)
    
    #-RED-KILONOVAE-----------------------------------------------
    colorGWr,deltamagGWr = Calculate_ColorDelta(dataGWr['g'],dataGWr['i'],delta_t,color_t)
    
    
    #-RISING COOLING ENVELOPE (16gkg)------------------------------

    colorIIb_e,deltamagIIb_e = Calculate_ColorDelta(dataIIb['g'],dataIIb['i'],delta_t,color_t)
    
    #-DECLINING COOLING ENVELOPE (16gkg)------------------------------
    colorIIb_l,deltamagIIb_l = Calculate_ColorDelta(dataIIb_l['g'],dataIIb_l['i'],delta_t,color_t)
    
    
    #-EARLY TYPE IA SN--(< 6 days from explosion)------------------
    #-SN2017cbv----------------------------------------------------
    colorIa_e,deltamagIa_e = Calculate_ColorDelta(dataIa['g'],dataIa['i'],delta_t,color_t)

    #-SN2011fe------------------------------------------------------
    colorIa2,deltamagIa2 = Calculate_ColorDelta(dataIa2['g'],dataIa2['i'],delta_t,color_t)
    
    
    #-NORMAL TYPE IA SN (> 10 days post-explosion-----------------_
    #-SN2017cbv----------------------------------------------------
    colorIa_l,deltamagIa_l = Calculate_ColorDelta(dataIa_l['g'],dataIa_l['i'],delta_t,color_t)
    
    #-SN2011fe------------------------------------------------------
    colorIa2_l,deltamagIa2_l = Calculate_ColorDelta(dataIa2_l['g'],dataIa2_l['i'],delta_t,color_t)
    
    
    #-FBOTS: RISING PHASE-----------------------------------------
    #-PS10bjp-----------------------------------------------------
    colorFB2,deltamagFB2 = Calculate_ColorDelta(dataFB2['g'],dataFB2['i'],delta_t,color_t)
    
    #-PS10ah-----------------------------------------------------
    colorFB3,deltamagFB3 = Calculate_ColorDelta(dataFB3['g'],dataFB3['i'],delta_t,color_t)
    
    #-PS12brf-----------------------------------------------------
    colorFB4,deltamagFB4 = Calculate_ColorDelta(dataFB4['g'],dataFB4['i'],delta_t,color_t)
    
    
    #-FBOTS: DECLINING PHASE--------------------------------------
    #-PS10ah------------------------------------------------------
    colorFB3b,deltamagFB3b = Calculate_ColorDelta(dataFB3b['g'],dataFB3b['i'],delta_t,color_t)
    
    #-PS12brf------------------------------------------------------
    colorFB4b,deltamagFB4b = Calculate_ColorDelta(dataFB4b['g'],dataFB4b['i'],delta_t,color_t)
    
    #-PS11qr--------------------------------------------------------
    colorFB6,deltamagFB6 = Calculate_ColorDelta(dataFB6['g'],dataFB6['i'],delta_t,color_t)
    
    #-AT2018cow------------------------------------------------------
    colorFB7,deltamagFB7 = Calculate_ColorDelta(dataFB7['g'],dataFB7['i'],delta_t,color_t)

    
    #-RAPIDLY DECLINING TYPE I------------------------------------
    #-SN2005ek----------------------------------------------------
    colorFB5,deltamagFB5 = Calculate_ColorDelta(dataFB5['g'],dataFB5['i'],delta_t,color_t)

    normals = []
    for i,sn in enumerate(sne):
        sn = sn[(sn['time-rel'] < 10*20) * (sn['time-rel'] > -7*20)]
        normals.append(Calculate_ColorDelta(sn['g'], sn['i'], delta_t, color_t))

    returnvalue =  [colorIa_e,deltamagIa_e,
                    colorIa2,deltamagIa2 ,
                    colorIa_l,deltamagIa_l,
                    colorIa2_l,deltamagIa2_l,
                    colorGW,deltamagGW ,
                    colorGWr,deltamagGWr,
                    colorIIb_e,deltamagIIb_e,
                    colorIIb_l,deltamagIIb_l,
                    colorFB2,deltamagFB2,
                    colorFB3,deltamagFB3,
                    colorFB4,deltamagFB4,
                    colorFB3b,deltamagFB3b,
                    colorFB4b,deltamagFB4b,
                    colorFB6,deltamagFB6,
                    colorFB7,deltamagFB7,
                    colorFB5,deltamagFB5, normals, len(normals)]

    return returnvalue
    



'''
        
tmp = readdata()

tmp2 = getcdt(data[0], data[1], data[2], data[3], data[4],
           delta_t=3.0, color_t=3.0)

plot(tmp2)

'''
from sklearn.svm import SVC
import pylab as pl

from sklearn.metrics import accuracy_score


dt1 = [1.5]#[0.5, 1.5, 3.5]#[0.5,1,1.5,2.5,3.5,4.5,5.5,6.5]
dt2 = [0.5]#[0,0.5,1]
data = readdata()

newdata = []
for i,d in enumerate(data[:-2]):
    newdata.append(d[d['time-rel'] < 15*24])

newdata.append(data[-2])
newdata.append(data[-1])               
               
    #d = d[d['time-rel'] < 30]
#pl.ion()

if __name__ == '__main__':
    print("here")
    kernel = 1.0 * RBF([1.0, 1.0])
    clf = GaussianProcessClassifier(kernel)
 
    for i in [(t1, t2) for t1 in dt1 for t2 in dt2]:
        print(i)
        tmp2 = getcdt(newdata,
                  delta_t=i[0], color_t=i[1])

        
        color = np.concatenate([tmp2[i*2] for i in range(15)])
        shape = np.concatenate([tmp2[i*2+1] for i in range(15)])
        #nsnIa = np.array([len(tmp2[i-2]) for i in range(4)]).sum()
        iacolor = np.concatenate([sn[0] for sn in tmp2[-2]])
        iashape = np.concatenate([sn[1] for sn in tmp2[-2]])
        N = len(color)
        randIa = np.random.randint(0, len(iacolor), N)
        plotcolor = ['IndianRed'] * N + ['k'] * N
        label = np.concatenate([np.zeros(N), np.ones(N)])
        indx = np.random.randint(0,len(label),size=int(N/2))

        print(len(color), N + len(iacolor), len(plotcolor), len(label))
        
        color = np.concatenate([color, iacolor[randIa]])
        shape = np.concatenate([shape, iashape[randIa]])

        X = np.array([shape, color]).T[indx]
        y = label[indx]
        
        xx = np.linspace(X.T[0].min(), X.T[0].max(), 100)
        yy = np.linspace(X.T[1].min(), X.T[1].max(), 100).T
        xx, yy = np.meshgrid(xx, yy)
        Xfull = np.c_[xx.ravel(), yy.ravel()]

        pl.plot(shape[-N:], color[-N:], '.', c='k', alpha=0.01)    
        pl.plot(shape[:-N], color[:-N], '.', c='IndianRed')
        pl.draw()
        #nsnIa = np.array([len(tmp2[i-2]) for i in range(4)]).sum()
        
        #label = np.zeros(len(color))
        #label[:nsnIa] = 1

        clf.fit(X, y)

        y_pred = clf.predict(X)
        accuracy = accuracy_score(y, y_pred)
        print("Accuracy (train) for %0.1f%% " % (accuracy * 100))

        # View probabilities:
        probas = clf.predict_proba(Xfull)
        n_classes = np.unique(y_pred).size
        
        imshow_handle = pl.imshow(probas[:, 0].reshape((100, 100)),
                                   extent=(X.T[0].min(),
                                            X.T[0].max(), X.T[1].min(),
                                            X.T[1].max()), 
                                           aspect='auto', origin='lower')
        idx = (y_pred == 0)
        if idx.any():
            pl.scatter(X[idx, 0], X[idx, 1], marker='o', c='w', edgecolor='k')

        pl.title("Probability")
        pl.show()
        input()
        
        
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
    print("Accuracy (train) for  %0.1f%% " %( accuracy * 100))

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
