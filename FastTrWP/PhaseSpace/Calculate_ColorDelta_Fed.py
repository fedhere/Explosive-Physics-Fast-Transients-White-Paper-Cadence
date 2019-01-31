'''Function to calculate the vectors for the delta_mag and color give a set of inputs.'''
import numpy as np

def Calculate_ColorDelta(g_array, i_array, delta_t, color_t):
    '''
    :param g_array: array of g-band magnitudes, already interpolated to 30 min timesteps
    :param i_array: array of i-band magnitudes, already interpolated to 30 min timesteps
    :param delta_t: change in time between subsequnt g-band measurements. Units of 0.5 hours.
    :param color_t: change in time between g and i band measurements for color. Units of 0.5 hours.
    :return: color, and deltamag arrays.
    '''
    #print(len(g_array))
    mag_step = np.int(delta_t / 0.5)
    color_step = np.int(color_t / 0.5)

    nmts = len(g_array) - mag_step
    ncts = len(g_array) - color_step    
    # color = list()
    #deltamag = list()

    recolor = np.zeros(nmts)
    redeltamag = np.zeros(ncts)

    #print(color_step, mag_step)
    #print(g_array[:10], i_array[:10])
    #for i in range(0, len(g_array) - mag_step):
    #    color.append(g_array[i] - i_array[i + color_step])
    #    deltamag.append(g_array[i] - g_array[i + mag_step])
    
    #color = np.array(color)
    #deltamag = np.array(deltamag)
    if color_step > 0.:
        recolor = g_array[:-color_step] -  i_array[color_step:]
    else:
        recolor = g_array -  i_array
    redeltamag = g_array[:-mag_step] - g_array[mag_step:]

    # pick only points with color and shape info (depends in the gaps)
    nmax = min(recolor.shape[0], redeltamag.shape[0]) 
    recolor = recolor[:nmax]
    redeltamag = redeltamag[:nmax]    

    #print(recolor.shape, redeltamag.shape)
    #print (color[:10], recolor[:10])
    #print("\n\n")

    #print (deltamag[:10], redeltamag[:10])
    
    return recolor, redeltamag
