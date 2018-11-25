# # KASEN MODEL FOR SN Ia COMPANION INTERACTION 

# ### this is code that takes a set of simulated spectra produced by Dan Kasen according to Kasen 2010 http://adsabs.harvard.edu/abs/2010ApJ...708.1025K and integrates them through the Kepler R filter to generate lightcurves of the first 10 days os a SN explosion.  $$ $$  The simulations reproduce the effect of a companion in a single degenerate SN Ia explosion. The companion is assumed to be in Roche lobe overflow, thus its distance is determined by its mass. Three companion cases are modelled:  $$ $$ a 2 $M_\odot$ Main Sequence (MS) star,  $$ $$ a 6 $M_\odot$ MS star, and  $$ $$ a 1 $M_\odot$ Red Giant (RG).  $$ $$  Note the important note about the lightcurve shape in the last paragraph of this notebook

# Time is sampled between 0.05 and 0.95 days after explosion with 0.1 days sampling intervals

import os,sys,glob
import numpy as np
import pylab as pl
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib.ticker import ScalarFormatter
formatter=ScalarFormatter()
formatter.set_scientific(True)
formatter.set_powerlimits((-2,2))





READONLY = True
#False
if READONLY:

        for f in ['r', 'g', 'i']:
                for dist in ['a5e11','a2e12','a2e13']:
                        ax = pl.figure().add_subplot(111)                
                        pd.read_csv("%s_%s_smoothed_baselineremoved.csv"%(dist,
                                        f)).plot(x="phase",
                                                 y="12.8386", ax=ax)
                        pd.read_csv("%s_%s_smoothed_baselineremoved.csv"%(dist,
                                        f)).plot(x="phase",
                                                 y="167.1616", ax=ax)
                        ax.set_title("%s %s"%(dist, f))
                        pl.show()
        sys.exit()
PLOTALL = False #True
PLOTS = True
# ### 40 viewing angles are simulated, between 0 and 180 degrees. $$ $$ 0 degree has the companion between us ant the WD. $$ $$ 180 has the WD between us and the companion. $$ $$ The viewing angle sampled regularly in cosine space, to give uniform elements of surface area on the sphere. 

angles = [167.1616, 157.6685, 151.0451, 145.5886,
          140.8052, 136.469, 132.4543, 128.6823,
          125.0997, 121.6683, 118.3595, 115.1508, 
          112.0244, 108.9657, 105.9621, 103.0030, 
          100.0787,  97.1808,  94.3013,  91.4326,  
          88.5675,  85.6989,  82.8193,  79.9214,  
          76.9972,  74.0380,  71.0345,  67.9757,  
          64.8494,  61.6407,  58.3318,  54.9004,  
          51.3179,  47.5459,  43.5312,  39.1950,  
          34.4115,  28.9550,  22.3317,  12.8386]
distcolor = {'a5e11': 'k', 'a2e12': 'r', 'a2e13': 'c'}
lastmu = angles[-2]
labels = {"a5e11": r"2 $M_\odot$ MS",
          "a2e12": r"6 $M_\odot$ MS",
          "a2e13": r"1 $M_\odot$ RG"}

# smoothing function, hanning window

def smooth(x,window_len=11):
        if x.ndim != 1:
            print("smooth only accepts 1 dimension arrays.")
            raise ValueError
        if x.size < window_len:
            print("Input vector needs to be bigger than window size.")
            raise ValueError
                
        if window_len<3:
                return x

        s = np.r_[2 * x[0] - x[window_len - 1::-1], x,
                  2 * x[-1] - x[-1:-window_len:-1]]

        w = eval('np.hanning(window_len)')
        y = np.convolve(w / w.sum(), s, mode='same')
        return y[window_len:-window_len+1]

if PLOTS:
        fig3 = {}
        fig4 = {}
        fig5 = {}        

# alphas and colors for plotting
dists = ['a5e11', 'a2e12', 'a2e13']
alphas = {}
for mu in angles:
    alphas[mu] = 1.0 - mu / max(angles)
allcolors=[ "Teal",  "YellowGreen", "Plum",  "SlateBlue", 
           "RoyalBlue", "Tomato", "RosyBrown", "Brown",  
           "SeaGreen",  "Crimson", "Cyan", "DarkBlue", 
           "DarkCyan", "DarkGreen", "DarkGray", "DarkGoldenRod", 
           "DarkKhaki", "DarkMagenta", "DarkOliveGreen", "DarkOrange", 
           "DarkOrchid", "DarkRed", "Olive","DarkSlateBlue", 
           "DarkSlateGray", "DarkTurquoise",  "Purple", "DeepPink", 
           "DeepSkyBlue", "DimGray", "DodgerBlue", "FireBrick", 
           "Turquoise", "ForestGreen", "Fuchsia", "OliveDrab", 
           "Gold", "GoldenRod", "Gray", "Green", "GreenYellow",
           "IndianRed", "Indigo",  "SteelBlue", "Khaki", "Lavender",  
           "LawnGreen","SpringGreen",  "Lime", "LimeGreen", "Magenta", 
           "Maroon", "MediumAquaMarine", "MediumBlue", "MediumOrchid", 
           "MediumPurple", "MediumSeaGreen", "MediumSlateBlue", 
           "MediumSpringGreen", "MediumTurquoise", "MediumVioletRed", 
           "MidnightBlue",  "Navy",  "OliveDrab", "Orange", "OrangeRed", 
           "Orchid",  "Peru", "Pink", "SlateBlue", "Red", "RosyBrown", 
           "RoyalBlue", "SaddleBrown", "Salmon", "SandyBrown", "Sienna", 
           "Silver", "SkyBlue", "Purple", "SlateGray", "SpringGreen", 
           "SteelBlue", "Tan", "Turquoise"]

myfilter = {}
filtfunc = {}
smoothed = {}
aspectra = {'w': np.arange(51, 14852, 100)}
lcvs = {'epoch': np.arange(0.05,10,0.1)}
smoothed = {}

allspectrafiles = {}
for dist in dists:
    allspectrafiles[dist] = {}
    spdir = 'companions/companion_' + dist
    for t in np.arange(0.05,10,0.1):
            allspectrafiles[dist][t] = np.array(glob.glob(spdir +
                                "/optical_t%06.2f_I3.spec.*"%t))

#loop over filters
for f in ['g', 'r', 'i']:
    if PLOTS:
        fig3[f] = pl.figure(figsize=(10,10))
        fig4[f] = pl.figure(figsize=(10,10))
        fig5[f] = pl.figure(figsize=(10,10))        
        
    filterfile = '%s.txt'%f 
    
    #used B filter for testing and comparing with Kasen 2010
    #filterfile = '/Users/fbianco/science/Dropbox/idl2/B_filter.dat'

    myfilter[f] = np.loadtxt(filterfile,unpack=True)
    # if filter files in nm convert to AA
    myfilter[f][0] *= 10
    filtfunc[f] = interp1d(myfilter[f][0], myfilter[f][1],
                         bounds_error=False, fill_value=0.0)
    if PLOTALL:
            fig0 = pl.figure()
            pl.plot(myfilter[f][0], myfilter[f][1], 'k--', label="%s"%f)
            pl.xlabel(r"wavelength (\AA)")
            pl.ylabel("transmission fraction")
            pl.legend(prop={'size':15})

    #pl.show()
    #continue
    # A distance string goes into the file names and identifies the companion.
    # It can be: 'a5e11', that is your 2Msun MS star, 'a2e12', your 6Msun MS star,
    # or 'a2e13', your Red Giant.
    lcvs[f] = {}

    #for each filter loop over companions
    for dist in dists:
        lcvs[f][dist] = {}
        smoothed[dist] = {}
        
        # #### filtering, integrating, plotting the spectra here

        print ("%s"%labels[dist])
        print ("each phase is plotted in a box in each plot arrays")
        print ("for each phase the spectrum at the narrowest " +
               "and largest viewing angles are plotted")
        print ("(only the narrowest is visible thought due to scaling)") 


        for mu in angles:
                lcvs[f][dist][mu] = []
                #print(mu)

        if PLOTALL:
                #fig1 unfiltered spectra
                fig1 = pl.figure(figsize=(15,15))
                fig2 = pl.figure(figsize=(15,15))          
                fig1.suptitle(r"%s unfiltered spectra %s"%(labels[dist], f), fontsize=20)
                fig2.suptitle(r"%s filtered spectra %s"%(labels[dist], f), fontsize=20)
                

        for t,epoch in enumerate(lcvs['epoch']):
                # for each epoch read in a spectrum for each angle
                for mu in angles:
                    thisf = [ff.strip()  for ff in allspectrafiles[dist][epoch] if
                             '%.4f'%mu in ff]
                    aspectra[mu] = np.loadtxt(thisf[0], usecols=(1,))
                if PLOTALL:
                        ax1 = fig1.add_subplot(10, 10, t+1)
                        ax2 = fig2.add_subplot(10, 10, t+1)            
                        ax1.text(300, 1500, "day: %.2f"%epoch)
                        ax2.text(300, 400, "day: %.2f"%epoch)
        
                for k in angles:
                        if PLOTALL:
                                if k == 'w':
                                        continue
                                if k == 167.1616:    
                                        ax1.plot(aspectra['w'] * 0.1,
                                                 aspectra[k], alpha=0.5,
                                                 color=allcolors[t%10])
                                        ax2.plot(aspectra['w'] * 0.1,
                                                 aspectra[k] * filtfunc[f](aspectra['w']),
                                                 alpha=0.5, color=allcolors[t%10])

                                if k == 12.8386:    
                                        ax1.plot(aspectra['w'] * 0.1,
                                                 aspectra[k],alpha=1.0,
                                                 color=allcolors[t%10])
                                        ax2.plot(aspectra['w'] * 0.1,
                                                 aspectra[k] * filtfunc[f](aspectra['w']),
                                                 alpha=1, color=allcolors[t%10])
                                
                        #integrating under the filter to get the lightcurve datapoint
                        lcvs[f][dist][k].append(sum(aspectra[k] *
                                        filtfunc[f](aspectra['w'])))
                if PLOTALL:
                        ax1.set_ylim(0, 2000)
                        ax1.set_yticks([0, 1000])
                        ax1.set_xticks([500, 1000])
                        #ax1.xaxis.set_major_formatter(formatter)
                        ax1.yaxis.set_major_formatter(formatter)
                        
                        ax2.plot(aspectra['w'] * 0.1, 500 * filtfunc[f](aspectra['w']),
                                 'k--', alpha=1)
                        #ax2.set_ylim(0,500)
                        ax2.set_yticks([0, 250])
                        ax2.set_xticks([500, 1000])
                        #ax2.xaxis.set_major_formatter(formatter)
                        ax2.yaxis.set_major_formatter(formatter)
                        #ax2.set_xlim([300,1000])

                        if (t) % 10 == 0:
                                ax1.set_ylabel("flux")
                                ax2.set_ylabel("flux")

                        if (t) / 10 == 9:
                                ax1.set_xlabel(r"wavelength ($\mu$)")
                                ax2.set_xlabel(r"wavelength ($\mu$)")

        if PLOTS:
                axfig3 = fig3[f].add_subplot(111)
                for i,mu in enumerate(angles):
                        axfig3.plot(lcvs['epoch'], lcvs[f][dist][mu], 
                            '%s'%distcolor[dist], alpha=(1.0-mu/max(angles)) / 2)
                axfig3.set_title(r"companion: %s star %s"%(labels[dist], f),
                                 fontsize=20)
                axfig3.set_ylabel("flux (arbitrary units?)", fontsize=15)
                axfig3.set_xlabel("days since explosion", fontsize=15);
                fig3[f].show()
        # printing file out

        theselcvs = {'phase': lcvs['epoch']}

        for mu in angles:
                theselcvs[mu] = lcvs[f][dist][mu]
        theselcvs = pd.DataFrame(theselcvs)

        theselcvs.to_csv("%s_%s.csv"%(dist, f), index=False)

        # ### now, the spectra are hella noisy, and so are the lightcurves!
        #  so i am smoothing them. smoothing window is ~1 day

        if PLOTS:
                axfig4 = fig4[f].add_subplot(111)                        

        for i,mu in enumerate(angles):
                smoothed[dist][mu] = smooth(np.array(lcvs[f][dist][mu]),
                                               window_len=11)

                if PLOTS:
                        axfig4.plot(np.array(lcvs['epoch']), smoothed[dist][mu],
                                    '%s--'%distcolor[dist], 
                                    alpha=1.0 - mu / max(angles))
        if PLOTS:
                axfig4.plot(np.array(lcvs['epoch']),
                            smoothed[dist][lastmu], '%s-'%distcolor[dist], 
                            alpha=1.0 - lastmu / max(angles),
                            label=labels[dist])

                pl.draw()
                pl.grid()
                axfig4.set_yscale('log')
                axfig4.set_ylim(-10, 3.0e4)
                pl.tick_params(axis='both', which='major', labelsize=15)
                pl.legend(loc=4, prop={'size':15})
                axfig4.set_ylabel("flux (arbitrary units?)", fontsize=15)
                axfig4.set_xlabel("days since explosion", fontsize=15)
                axfig4.set_title("smoothed lightcurves (log scale) - %s"%f, fontsize=20)



        #print smoothed out
        theselcvssmooth = {'phase': lcvs['epoch']}
        for mu in angles:
                theselcvssmooth[mu] = smoothed[dist][mu]
        theselcvssmooth = pd.DataFrame(theselcvssmooth)
        theselcvssmooth.to_csv("%s_%s_smoothed.csv"%(dist, f), index=False)
    if PLOTS:
            axfig5 = fig5[f].add_subplot(111)
    for dist in dists:
        for i,mu in enumerate(angles[:-1]):

                if PLOTS:
                        axfig5.plot(np.array(lcvs['epoch']), 
                                smoothed[dist][mu] / smoothed['a5e11'][angles[0]],
                                '%s-'%distcolor[dist], 
                                    alpha=(1.0-mu/(angles[0])+0.01)/2)
        if PLOTS:
                axfig5.plot(np.array(lcvs['epoch']), 
                    smoothed[dist][lastmu] / smoothed['a5e11'][angles[0]],
                    '%s-'%distcolor[dist], 
                    alpha=1.0-lastmu/(angles[0])+0.01, label=labels[dist])

                pl.draw()
                pl.grid()
                #pl.yscale('log')
                pl.ylim(-0.3,60)
                pl.tick_params(axis='both', which='major', labelsize=15)
                #pl.legend(prop={'size':15})
                pl.ylabel("relative flux",fontsize=15)
                pl.xlabel("days since explosion",fontsize=15)
                pl.title("excess only time series %s"%f,fontsize=20)

        theselcvsnormed = {'phase': lcvs['epoch']}
        for mu in angles:
                theselcvsnormed[mu] = smoothed[dist][mu] / smoothed['a5e11'][angles[0]]
        theselcvsnormed = pd.DataFrame(theselcvsnormed)

        theselcvsnormed.to_csv("%s_%s_smoothed_baselineremoved.csv"%(dist, f),
                               index=False)
pl.show()

