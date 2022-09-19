import numpy as np
import matplotlib.pyplot as plt
from Lu_decomp import LUdecompose
from scipy.interpolate import CubicSpline

# This file answers question 3 from the assignment

def lagrange(x,y):

    """""
    :param x: Input array for the x-axis
    :param y: Input array for the y-axis
    :return array_x: outputs an array of x values using lagrange polynomial interpolation
    :return array_y: outputs an array of y values using lagrange polynomial interpolation

    This function performs lagrange interpolation on our data set x and y.
    """

    x_len = len(x)                                    # Obtain the length of x
    array_x = np.linspace(x[0], x[x_len - 1], 1000)   # create array of x values using our original x, with 1000 samples
    array_y = []                                      # Empty list for the y values
    for k in range(len(array_x)):                     # iterate through the x array values
        val = array_x[k]                              # specific value in the array_x
        tot = 0
        for i in range(x_len):                        # Implementation of lagrange polynomial function
            store = 1
            for j in range(x_len):
                if i != j:                            # for different values
                    store = store * ((val - x[j])/ (x[i]-x[j]))         # Lagrange equation
            store = store * y[i]                      # update store
            tot = tot + store                         # calculating y value
        array_y.append(tot)                           # place the y value into the y array
    return array_x, array_y


def spline(x,y, acc):

    """"
    :param x: input array of x values
    :param y:  input array of y values
    :param acc: The number of samples which also directs the accuracy of the spline

    :return t:  Returns an array of time values
    :return s: Returns the spline values, which is essentially the y array for plotting

    This function performs a cubic spline on a data set x and y using the natural spline boundary. The function
    LUdecompose is called to solve the martrix, which is in the file Lu_decomps
    """

    x_dim = np.shape(x)[0]
    t = np.linspace(x[0], x[x_dim-1], acc)  # create a matrix for the time values
    A = np.zeros((x_dim, x_dim))       # create a matrix of the values that are used to calculate the second derivatives
    b = np.zeros((x_dim, 1))           # create a matrix for the b values
    s = np.zeros((acc, 1))             # array to store the splines
    for i in range(x_dim):             # iterate through the number of dimensions
        if i == 0:
            A[0][0] = 1                          # setting the first term in A to 1
            b[0] = 0                             # setting first term in b to 0
        elif i == x_dim - 1:
            A[x_dim - 1][x_dim - 1] = 1          # setting last term in A to 1
            b[x_dim - 1] = 0                     # setting the last term in b to 0
        else:
            A[i][i - 1] = (x[i] - x[i - 1]) / 6  # Adding the values for the matrix
            A[i][i] = (x[i + 1] - x[i - 1]) / 3
            A[i][i + 1] = (x[i + 1] - x[i]) / 6
            b[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])  # adding b values

    f = LUdecompose(A, b)             # Using my LU decomposition,found the matrix values for my second derivatives
    for i in range(len(t)):  # iterate over all values in t
        for j in range(x_dim-1):  # iterate over the values in x
            if x[j] <= t[i] <= x[j + 1]:  # calculate the spline for this region
                a = (x[j + 1] - t[i]) / (x[j + 1] - x[j])  # calculating the coefficients
                b = (t[i] - x[j]) / (x[j + 1] - x[j])
                c = (1 / 6) * (a ** 3 - a) * (x[j + 1] - x[j]) ** 2
                d = (1 / 6) * (b ** 3 - b) * (x[j + 1] - x[j]) ** 2
                s[i] = a * y[j] + b * y[j + 1] + c * f[j] + d * f[j + 1]
    return t , s


# Define the variables that will be used in the function
x_input = np.array([-0.75, -0.5, -0.35, -0.1, 0.05, 0.1, 0.23, 0.29, 0.48, 0.6, 0.92, 1.05, 1.5])
y_input = np.array([0.1, 0.3, 0.47, 0.66, 0.60, 0.54, 0.3, 0.15, -0.32, -0.54, -0.6, -0.47, -0.08])

# The following code runs the function lagrange and cubic spline function and plots the function
lag_array_x, lag_array_y = lagrange(x_input , y_input)
cubic_x_array,cubic_spline_array = spline(x_input, y_input, 10000)
plt.figure(1)
plt.xlabel("X Axis")
plt.ylabel("Y Axis")
plt.plot(lag_array_x, lag_array_y, label = 'Lagrange Polynomial ')
plt.plot(cubic_x_array, cubic_spline_array, label = ' Programmed Cubic Spline ')
plt.plot(x_input, y_input, label = 'Original Data', ls = ' ',  marker = 'x')
plt.legend()
plt.show()


# The following code runs the function spline and plots the function
# cubic_x_array,cubic_spline_array = spline(x_input, y_input, 10000)
#
# plt.figure(2)
# plt.xlabel("X Axis")
# plt.ylabel("Y Axis")
# plt.plot(cubic_x_array, cubic_spline_array, label = ' Programmed Cubic Spline ')
# plt.plot(x_input, y_input, label = 'Original Data', marker = 'x')
# plt.legend()
# plt.show()

# Following code plots the function of an in-built cubic spline to compare to my own function
c = CubicSpline(x_input,  y_input, bc_type = 'natural')      # using a cubic spline function from Scipy

plt.figure(2)
plt.xlabel("X Axis")
plt.ylabel("Y Axis")
plt.plot(cubic_x_array, cubic_spline_array, label = ' Programmed Cubic Spline ', linewidth = 7, color = 'c')
plt.plot(cubic_x_array, c(cubic_x_array) , label = 'Scipy Cubic Spline', color = 'm', linewidth = 3)
plt.plot(x_input, y_input, label = 'Original Data', marker = 'x', color = 'orange')
plt.legend()
plt.show()


#since the cubioc spline works then it should be fine

