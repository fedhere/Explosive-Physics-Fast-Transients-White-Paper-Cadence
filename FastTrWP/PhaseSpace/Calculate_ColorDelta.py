'''Function to calculate the vectors for the delta_mag and color give a set of inputs.'''
import numpy as np

def Calculate_ColorDelta(g_array,i_array,delta_t,color_t):
    '''
    :param g_array: array of g-band magnitudes, already interpolated to 30 min timesteps
    :param i_array: array of i-band magnitudes, already interpolated to 30 min timesteps
    :param delta_t: change in time between subsequnt g-band measurements. Units of 0.5 hours.
    :param color_t: change in time between g and i band measurements for color. Units of 0.5 hours.
    :return: color, and deltamag arrays.
    '''

    step = np.int(delta_t / 0.5)
    color_step = np.int(color_t / 0.5)

    color = list()
    deltamag = list()

    for i in range(0, len(g_array) - step):
        color.append(g_array[i] - i_array[i + color_step])
        deltamag.append(g_array[i] - g_array[i + step])

    color = np.array(color)
    deltamag = np.array(deltamag)

    return color, deltamag