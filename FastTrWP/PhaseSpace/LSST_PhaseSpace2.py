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
from bokeh.models import Range1d, Label

from Calculate_ColorDelta import Calculate_ColorDelta
import glob
import sys
#Bokeh setup
from bokeh.plotting import figure,output_file,show, ColumnDataSource, gridplot

#Set some initial parameters:
delta_t = 4.0 #set in terms of number of hours. Should always be in units of 0.5 hours.
color_t =  0.5 #set in hours. Must either be 0 or 0.5 (for 30 minutes)

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

p = figure(title="∆T1 = 0.5 hours; ∆T2 = 4 hours", x_axis_label='change in g-band magnitude (mag)', y_axis_label='(g - i) color (mag)',tools=tools)


p.patch([0.04,0.62,0.62,0.04],[-1,-1,1,1],alpha=0.1,fill_color='MediumSpringGreen')
p.patch([-0.03,-0.4,-0.4,-0.03],[-1,-1,0.2,0.2],alpha=0.1,fill_color='LightSeaGreen')
p.patch([-0.03,-0.4,-0.4,-0.03],[0.3,0.3,2.5,2.5],alpha=0.1,fill_color='Red')


#Need to bulk read in Fed's Type Ia here so that they can be plotted immediately and don't separately save all of them.
paths = glob.glob('./IAs/snIa*.csv')
print (paths)
for path in paths:
    data = ascii.read(path)
    index = [(data['time-rel'] > -5*24) & (data['time-rel'] < 20*24)]
    gIa = data['g'][index]
    iIa = data['i'][index]
    colorIa,deltamagIa = Calculate_ColorDelta(gIa,iIa,delta_t,color_t)
    p.circle(deltamagIa[0::12],colorIa[0::12],legend='Type Ia (normal)',size=4.5,fill_color='Black',fill_alpha=0.02,line_alpha=0)
    #only plot 1/day?

#Need to bulk read in Fed's Type Ia here so that they can be plotted immediately and don't separately save all of them.
#paths = glob.glob('./LightCurves/Formatted_others/sn*.dat')
#for path in paths:
#    data = ascii.read(path)
#    colorIa,deltamagIa = Calculate_ColorDelta(data['g'],data['i'],delta_t,color_t)
#    p.square(deltamagIa,colorIa,legend='Type Ia',line_width=0.5,size=2.5,fill_color='Orange',line_color='Black')


#Type Ia SN: Normal
#p.line(deltamagIa_l,colorIa_l,legend='Type Ia (peak)',line_width=0.5,line_color='Black')
#p.circle(deltamagIa_l,colorIa_l,legend='Type Ia (peak)',line_width=0.5,size=2.5,fill_color='Grey',line_color='Black')
#p.line(deltamagIa2_l,colorIa2_l,legend='Type Ia (peak)',line_width=0.5,line_color='Black')
p.square(deltamagIa2_l[0::12],colorIa2_l[0::12],legend='Type Ia (normal)',line_width=0.5,size=4.5,fill_color='Black',line_alpha=0,fill_alpha=0.01)


###########################################
#Kilonovae
#p.line(deltamagGW,colorGW,legend='Blue Kilonova',line_width=0.5,line_color='Blue')
p.triangle(deltamagGW[0:125],colorGW[0:125],legend='Blue Kilonova',line_width=0.5,size=8.5,fill_color='ForestGreen',line_color='Black')
#p.line(deltamagGWr,colorGWr,legend='Red Kilonova',line_width=0.5,line_color='Red')
p.triangle(deltamagGWr[0:125],colorGWr[0:125],legend='Red Kilonova',line_width=0.5,size=8.5,fill_color='IndianRed',line_color='Black')

p.triangle(deltamagGW[125:-400:12],colorGW[125:-400:12],legend='Blue Kilonova',line_width=0.5,size=8.5,fill_color='ForestGreen',line_color='Black')
#p.line(deltamagGWr,colorGWr,legend='Red Kilonova',line_width=0.5,line_color='Red')
p.triangle(deltamagGWr[125:-400:12],colorGWr[125:-400:12],legend='Red Kilonova',line_width=0.5,size=8.5,fill_color='IndianRed',line_color='Black')


#Rapidly declining Type I SN
#p.line(deltamagFB5[0::12],colorFB5[0::12],legend='Fast Decliner',line_width=0.5,line_color='Black')
p.circle(deltamagFB5[0::12],colorFB5[0::12],legend='Fast Decliner',line_width=0.5,size=8.5,fill_color='Orange',line_color='Black')


#FBOTS:
#p.line(deltamagFB2[0::6],colorFB2[0::6],legend='Fast Transient',line_width=0.5,line_color='Black')
p.circle(deltamagFB2[0::8],colorFB2[0::8],legend='Fast Blue Transient',line_width=0.5,size=9.5,fill_color='Cyan',line_color='Black')
#p.line(deltamagFB3[0::6],colorFB3[0::6],legend='Fast Transient',line_width=0.5,line_color='Black')
p.circle(deltamagFB3[0::8],colorFB3[0::8],legend='Fast Blue Transient',line_width=0.5,size=9.5,fill_color='Cyan',line_color='Black')
#p.line(deltamagFB4[0::6],colorFB4[0::6],legend='Fast Transient',line_width=0.5,line_color='Black')
p.circle(deltamagFB4[0::8],colorFB4[0::8],legend='Fast Blue Transient',line_width=0.5,size=9.5,fill_color='Cyan',line_color='Black',fill_alpha=0.5)


#p.line(deltamagFB3b[0::6],colorFB3b[0::6],legend='Fast Transient',line_width=0.5,line_color='Black')
p.circle(deltamagFB3b[0::12],colorFB3b[0::12],legend='Fast Blue Transient',line_width=0.5,size=9.5,fill_color='Cyan',line_color='Black',fill_alpha=0.5)
#p.line(deltamagFB4b[0::6],colorFB4b[0::6],legend='Fast Transient',line_width=0.5,line_color='Black')
p.circle(deltamagFB4b[0::12],colorFB4b[0::12],legend='Fast Blue Transient',line_width=0.5,size=9.5,fill_color='Cyan',line_color='Black',fill_alpha=0.5)

#p.line(deltamagFB6,colorFB6,legend='Fast Transient',line_width=0.5,line_color='Black')
p.circle(deltamagFB6[0::12],colorFB6[0::12],legend='Fast Blue Transient',line_width=0.5,size=10.5,fill_color='Cyan',line_color='Black',fill_alpha=0.5)
#p.line(deltamagFB7,colorFB7,legend='Fast Transient',line_width=0.5,line_color='Black')
p.circle(deltamagFB7[0::12],colorFB7[0::12],legend='Fast Blue Transient',line_width=0.5,size=10.5,fill_color='Cyan',line_color='Black',fill_alpha=0.5)


#Type IIb SN
#p.line(deltamagIIb_e,colorIIb_e,legend='Infant CCSN (Cooling Env.)',line_width=0.5,line_color='Black')
p.diamond(deltamagIIb_e,colorIIb_e,legend='Infant CCSN',line_width=0.5,size=10.5,fill_color='mediumblue',line_color='Black')

#p.line(deltamagIIb_l[0::6],colorIIb_l[0::6],legend='Infant CCSN (Cooling Env.)',line_width=0.5,line_color='Black')
p.diamond(deltamagIIb_l[0::6],colorIIb_l[0::6],legend='Infant CCSN',line_width=0.5,size=10.5,fill_color='mediumblue',line_color='Black')

#Type Ia SN: Early
#p.line(deltamagIa_e[0::6],colorIa_e[0::6],legend='Infant Type Ia',line_width=0.5,line_color='Black')
p.square(deltamagIa_e[0::8],colorIa_e[0::8],legend='Infant Type Ia',line_width=0.5,size=8.5,fill_color='Indigo',line_color='Black')
#p.line(deltamagIa2[0::6],colorIa2[0::6],legend='Infant Type Ia',line_width=0.5,line_color='Black')
p.square(deltamagIa2[0::12],colorIa2[0::12],legend='Infant Type Ia',line_width=0.5,size=8.5,fill_color='Indigo',line_color='Black')


p.legend.location = 'top_right'
p.xaxis.axis_label_text_font_size = '14pt'
p.xaxis.axis_label_text_font_style = 'normal'
p.yaxis.axis_label_text_font_size = '14pt'
p.yaxis.axis_label_text_font_style = 'normal'
p.yaxis.major_label_text_font_size = "12pt"
p.xaxis.major_label_text_font_size = "12pt"
p.x_range=Range1d(-0.46,0.62)
p.y_range=Range1d(-1.2,2.8)

p.title.align = "center"
p.title.text_font_size= "14pt"

#Labels:
blue = Label(x=-0.41,y=-0.6,text='BLUE',render_mode='canvas',angle=90,angle_units='deg',text_font_size='13pt',text_font_style='bold')
red = Label(x=-0.41,y=1.8,text='RED',render_mode='canvas',angle=90,angle_units='deg',text_font_size='13pt',text_font_style='bold')
decline = Label(x=-0.3,y=2.6,text='DECLINING',render_mode='canvas',text_font_size='12pt',text_font_style='bold')
rise = Label(x=0.1,y=2.6,text='RISING',render_mode='canvas',text_font_size='12pt',text_font_style='bold')
p.add_layout(blue)
p.add_layout(red)
p.add_layout(rise)
p.add_layout(decline)
p.line([-0.3,-0.1],[2.6,2.6],line_color='Black',line_width=4)
p.line([0.1,0.22],[2.6,2.6],line_color='Black',line_width=4)
p.line([-0.41,-0.41],[-0.6,-0.2],line_color='Black',line_width=4)
p.line([-0.41,-0.41],[1.8,2.1],line_color='Black',line_width=4)

fast_rise = Label(x=0.18,y=0.6,text='infant SN & fast risers',render_mode='canvas',text_font_size='10pt',text_font_style='italic')
p.add_layout(fast_rise)
fast_rise1 = Label(x=-0.37,y=0.0,text='cooling envelope &',render_mode='canvas',text_font_size='10pt',text_font_style='italic')
fast_rise2 = Label(x=-0.37,y=-0.2,text='blue fast decliners',render_mode='canvas',text_font_size='10pt',text_font_style='italic')
p.add_layout(fast_rise1)
p.add_layout(fast_rise2)
fast_rise1 = Label(x=-0.32,y=1.4,text='fading kilonova &',render_mode='canvas',text_font_size='10pt',text_font_style='italic')
fast_rise2 = Label(x=-0.32,y=1.2,text='red fast decliners',render_mode='canvas',text_font_size='10pt',text_font_style='italic')
p.add_layout(fast_rise1)
p.add_layout(fast_rise2)
#ACTUAL PLOT:

#Glyphs:
#Infant SN and Rapidly-Rising Transients
#Fading Cooling Envelope Emission and Fast Blue Transients
#Fading KN and Rapid red Transients

show(p)
