''' Create plot for delta mag and color for set of transients for science justification.

Use:
(1) Type Ia blue (17cbv)-done
(2) Type Ia red (11fe)-done
(3) KN Blue (GW170817)-done
(4) KN Red (Model of Gw170817 w/o blue)-done
(5) Type IIb cooling curve 16gkg done. Add more if time.
(6) FBOTs [gold templates and 18cow]. done
(8) Early type II light curves.
(9) Larger sample of 'normal' events. 

This currenly makes plots for g and i band colors, to match obssims that we have.'''


from astropy.io import ascii
import numpy as np
from Calculate_ColorDelta import Calculate_ColorDelta

#Bokeh setup
from bokeh.plotting import figure,output_file,show, ColumnDataSource, gridplot

#Set some initial parameters:
delta_t = 4.0 #set in terms of number of hours. Should always be in units of 0.5 hours.
color_t = 0. #set in hours. Must either be 0 or 0.5 (for 30 minutes)

step = np.int(delta_t/0.5)
color_step = np.int(color_t/0.5)

#----------------------------------------------------------------------
#-READ IN TEMPLATES----------------------------------------------------
#----------------------------------------------------------------------

#-KILONOVAE------------------------------------------------------------
dataGW = ascii.read('LightCurves/Formatted/GW170817.dat')
colorGW,deltamagGW = Calculate_ColorDelta(dataGW['g'],dataGW['i'],delta_t,color_t)

#-RED-KILONOVAE-----------------------------------------------
dataGWr = ascii.read('LightCurves/Formatted/GW170817red.dat')
colorGWr,deltamagGWr = Calculate_ColorDelta(dataGWr['g'],dataGWr['i'],delta_t,color_t)


#-RISING COOLING ENVELOPE (16gkg)------------------------------
dataIIb = ascii.read('LightCurves/Formatted/SN2016gkg_poly_e.dat')
colorIIb_e,deltamagIIb_e = Calculate_ColorDelta(dataIIb['g'],dataIIb['i'],delta_t,color_t)

#-DECLINING COOLING ENVELOPE (16gkg)------------------------------
dataIIb_l = ascii.read('LightCurves/Formatted/SN2016gkg_poly_l.dat')
colorIIb_l,deltamagIIb_l = Calculate_ColorDelta(dataIIb_l['g'],dataIIb_l['i'],delta_t,color_t)


#-EARLY TYPE IA SN--(< 6 days from explosion)------------------
#-SN2017cbv----------------------------------------------------
dataIa = ascii.read('LightCurves/Formatted/SN2017cbv_poly_e.dat')
colorIa_e,deltamagIa_e = Calculate_ColorDelta(dataIa['g'],dataIa['i'],delta_t,color_t)

#-SN2011fe------------------------------------------------------
dataIa2 = ascii.read('LightCurves/Formatted/SN2011fe_poly_e.dat')
colorIa2,deltamagIa2 = Calculate_ColorDelta(dataIa2['g'],dataIa2['i'],delta_t,color_t)


#-NORMAL TYPE IA SN (> 10 days post-explosion-----------------_
#-SN2017cbv----------------------------------------------------
dataIa_l = ascii.read('LightCurves/Formatted/SN2017cbv_poly_l.dat')
colorIa_l,deltamagIa_l = Calculate_ColorDelta(dataIa_l['g'],dataIa_l['i'],delta_t,color_t)

#-SN2011fe------------------------------------------------------
dataIa2_l = ascii.read('LightCurves/Formatted/SN2011fe_poly_l.dat')
colorIa2_l,deltamagIa2_l = Calculate_ColorDelta(dataIa2_l['g'],dataIa2_l['i'],delta_t,color_t)


#-FBOTS: RISING PHASE-----------------------------------------
#-PS10bjp-----------------------------------------------------
dataFB2 = ascii.read('LightCurves/Formatted/PS10bjp_rise.dat')
colorFB2,deltamagFB2 = Calculate_ColorDelta(dataFB2['g'],dataFB2['i'],delta_t,color_t)

#-PS10ah-----------------------------------------------------
dataFB3 = ascii.read('LightCurves/Formatted/PS10ah_rise.dat')
colorFB3,deltamagFB3 = Calculate_ColorDelta(dataFB3['g'],dataFB3['i'],delta_t,color_t)

#-PS12brf-----------------------------------------------------
dataFB4 = ascii.read('LightCurves/Formatted/PS12brf_rise.dat')
colorFB4,deltamagFB4 = Calculate_ColorDelta(dataFB4['g'],dataFB4['i'],delta_t,color_t)


#-FBOTS: DECLINING PHASE--------------------------------------
#-PS10ah------------------------------------------------------
dataFB3b = ascii.read('LightCurves/Formatted/PS10ah_decline.dat')
colorFB3b,deltamagFB3b = Calculate_ColorDelta(dataFB3b['g'],dataFB3b['i'],delta_t,color_t)

#-PS12brf------------------------------------------------------
dataFB4b = ascii.read('LightCurves/Formatted/PS12brf_decline.dat')
colorFB4b,deltamagFB4b = Calculate_ColorDelta(dataFB4b['g'],dataFB4b['i'],delta_t,color_t)

#-PS11qr--------------------------------------------------------
dataFB6 = ascii.read('LightCurves/Formatted/PS11qr_decline.dat')
colorFB6,deltamagFB6 = Calculate_ColorDelta(dataFB6['g'],dataFB6['i'],delta_t,color_t)

#-AT2018cow------------------------------------------------------
dataFB7 = ascii.read('LightCurves/Formatted/AT2018cow_decline.dat')
colorFB7,deltamagFB7 = Calculate_ColorDelta(dataFB7['g'],dataFB7['i'],delta_t,color_t)


#-RAPIDLY DECLINING TYPE I------------------------------------
#-SN2005ek----------------------------------------------------
dataFB5 = ascii.read('LightCurves/Formatted/SN2005ek_decline.dat')
colorFB5,deltamagFB5 = Calculate_ColorDelta(dataFB5['g'],dataFB5['i'],delta_t,color_t)



#------------------------------------------------------------------------------
#Set up Bokeh plot to interactively view the results:
#-----------------------------------------------------------------------------
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


show(p)