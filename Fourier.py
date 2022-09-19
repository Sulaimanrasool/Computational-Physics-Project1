import numpy as np
import math as mp
import matplotlib.pyplot as plt
import scipy as sp

# This file answers question 4 from the assignment

def convolve(t, y):
    """

    :param t: Value for time, which will act as the upper time limit (positive t) and lower time limit (negative t)
    :param y: Number of samples between -t to t

    :return h: The signal function, hard coded in convolve and padded
    :return g: The response function, hard coded in convolve and padded
    :return t: Array of time values for plotting
    :return con: The convolution of the signal and response function

    This function finds the convolution of the signal function with a response function
    """
    time = np.linspace(-t,t, y)                                 # create a time array

    h = np.zeros(y)                                          # array of zeros
    g = np.zeros(y)

    for i in range(y) :                                      # Creating the functions
        eq = (1/(mp.sqrt(2*np.pi))) * mp.exp((-time[i]**2)/4)   # making the response function
        g[i] = eq                                            # Appending the values
        if time[i] < 5 or time[i] > 7 :                   # Making the signal function (it is a step single square wave
            h[i] = 0
        else:
            h[i] = 4                                      # Value of the square wave

    np.pad(h, (0, y - 1), 'constant')                     # pad out h with zeros
    np.pad(g, (0, y - 1), 'constant')                     # pad out the g with zeros
    F1 = np.fft.fft(h, norm='ortho')                      # Fourier transform of signal function, normalise
    F2 = np.fft.fft(g, norm='ortho')                      # Fourier transform of response function, normalise
    con = np.fft.ifft(F2 * F1)                            # Inverse of F1 * F2, which is the convolution
    return h, g, time, con


# Initialise the variables
t = 20
number_of_samples = 1000

# Call the convolve function and obtain the variables for the plots
signal_function,responsive_function, time_array, convolution_array = convolve(t, number_of_samples)
convolve_shift = np.fft.fftshift(convolution_array).real  # Shifts the convolution

# Plot for the signal and response function
plt.xlabel("Time[s]")
plt.ylabel("Signal ")
plt.plot(time_array, signal_function, label = 'Signal Function' )
plt.plot(time_array, responsive_function, label = 'Response Function')
plt.legend()
plt.show()


# Plot for the Convolution
plt.xlabel("Time[s]")
plt.ylabel("Signal ")
plt.plot(time_array, convolve_shift, label = 'Convolution' )
plt.legend()
plt.show()


