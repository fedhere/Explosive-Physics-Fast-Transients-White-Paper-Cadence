# identify transients

MLG: identify literature and existing data of the types of features and LC that exhibit what we want to find, and use that data to define a diagnostic; timescale for the color changes too, and keeping in mind that you want less filter changes in 30 min. how fast do the light curves evolve, and is 30 minutes too short to see evolution. recheck fed's plots for the dm/dt, and whether intra night or next-night is better, and in what filter (same/different)

 - SN early bumps from shock breakout (CC), companions (Ia), stellar envelope and circumstellar material (IIb)
 - Kilonova
 - Microlensing?? 
 - .Ia?
 - Maria's objects


# define filters and time gaps and their tollerance

3 obs in 2 filters, 2 in different filters separated by <45 min, 1 in one of those filters separated by (?)>1.5 hour

and the next visit - do we need to come back in the same filter? 

extragalactic fields that are being included in rolling cadence

and maybe some GP fields (microlensing events on short timescales; binary stars and BH, not exoplanets)

TO DEFINE: 
filter pairs
FOR EACH FILTER PAIR:
    max filter gap 1 (for color) 
    min filter gap 2 (for shape) 
    max filter gap 2 (for timely ID) 
# metrics

describe and attach metrics for different science cases
The metric is the fraction of events for which the color and risetime and risetime is constrained (within some accuracy). 
Different science cases may have different input lightcurves and different gap requirements

first a diagnostic metric

then maybe a LC-detection metric



1 - given dt1 dt2 ranges and filter pairs counts the fraction of fields observed like that 

2 - given dt1 dt2 ranges and filter pairs  and event injection (monte carlo injection of template lightcurves) counts the fraction of detected events in all three observations within x days of the transient lifetime 

3  -  given dt1 dt2 ranges and filter pairs  and event injection (monte carlo injection of template lightcurves) for a couple of test cases (interacting SN Ia, KN, microlensing?) retrieves the lightcurve and measures the classification/characterization success

once we have 3 we can inject whatever we want and test wit whatever filters we want. so we know what tolerances in dt1 dt2 and filter pairs work for what science case. this can also be done with metric 2 if we assess what filter pairs and dt1-dt2 ranges are informative for a given model. The result of this can be summarized in a table in the technical details section (since we are running out of space in the sci just session)

# figures
we should also consider figures like figures 6.4, 6.5, 6.6 in here