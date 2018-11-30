import sys,os,glob, inspect
import numpy as np
from snclasses import mysn
import pyfits as pf

bands =['U','up','B','V','R','I','rp','ip', 'J','H','K']
if __name__=="__main__":

    if len(sys.argv)<2:
        print "Usage: python readsn_fitstable.py cfafitsfilename"
    snname=sys.argv[1].split('/')[-1].replace('fits','')
    f=pf.open(sys.argv[1])
    cols=f[1].columns
    allcols=cols.names
    tb=f[1].data
    thissn=mysn(tb.field('SNname')[0])
    thissn.type =tb.field('SNtype')[0]
    thissn.Vmax=tb.field('Vmaxdate')[0]
    try:
        float(thissn.Vmax)
        thissn.flagmissmax=False
    except:
        pass
    
    thissn.camsystem=tb.field('pipeversion')[0]
    for b in bands:
        if b in allcols:
            indx = np.where(tb.field(b+'epochs') > 0)[0]
            print b, indx
            thissn.filters[b[0]]=len(indx)
            thissn.photometry[b[0]]={}
            thissn.photometry[b[0]]['mjd']=np.array(tb.field(b+'epochs'))[indx]
            thissn.photometry[b[0]]['mag']=np.array(tb.field(b))[indx]
            thissn.photometry[b[0]]['dmag']=np.array(tb.field('d'+b))[indx]
            thissn.photometry[b[0]]['natmag']=np.array(tb.field(b+'_nat'))[indx]
            thissn.photometry[b[0]]['camsys']=np.array(tb.field(b+'pipeversion'))[indx]            
            
    thissn.printsn(photometry=True)
    thissn.plotsn(photometry=True,show=True, offsets=True)
