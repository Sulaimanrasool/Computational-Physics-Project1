import numpy as np
import matplotlib.pyplot as plt

# This file answers question 5 of the assignment

def volt_func(function, vi, r, c, n, t):

    """
    :param function: provide a string 'discharge' , 'square1', 'square2'. Specifies the type of fucntion to apply
    :param vi: Initial input voltage
    :param r:  Resistance
    :param c:  Capacitance
    :param n:  number of samples
    :param t:  maximum time limit

    :return voltage_array: returns an array of input voltages depending on the type of function
    :return tim: returns an array of time
    :return tow : returns the time constant

    This function returns a voltage array depending on the type of function being specified. The 'discharge' input sets
     the voltage in after t > 0 as 0. The 'square1' input creates an array of input voltages that act as an input
     square wave with the time period T1. The  'square2' input creates an array of input voltages that act as an input
     square wave with the time period T2.
    """
    tim = np.linspace(0, t, n )                          # create the array of time
    tow = r * c                                         # time constant
    v_d = []                                            # discharge function voltage array
    v_squ1 = []                                         # first square wave voltage array
    v_squ2 = []                                         # second square wave voltage array
    T1 = 2 * r * c                                      # Period definition for square1
    T2 = (r*c)/2                                        # Period definition for square2

    if function == 'discharge' :                        # creating array for the discharge
        for i in range(n):
            k = 0
            v_d.append(k)                               # append the values to the discharge voltage array
        return(v_d,tim, tow)

    if function == 'square1':                         # creating array for the square1 input
        for i in tim:
            if 0 <= i % T1 < T1/2 :                     # Allows for voltages to be calculated across the whole time span
                v_squ1.append(0)                       # Equation for trough
            if (T1/2) <= i % T1 < T1:                      # Equation for peak
                v_squ1.append(vi)
        return(v_squ1, tim, tow)

    if function == 'square2':                          # Creating an array for square2 input
        for i in tim:
            if 0 <= i % T2 < T2 / 2:                    # Inputs values for the trough into the list
                v_squ2.append(0)
            if (T2 / 2) <= i % T2 < T2:                    # Inputs values for the peak into the list
                v_squ2.append(vi)

        return (v_squ2, tim, tow)


def runge_kutta(t_array, vi, tow , n , vo):
    """
    :param t_array:   time array supplied created when obtaining voltage_in_array
    :param vi:  Array of voltage inputs depending on what input was supplied in volt_func
    :param tow: Time constant calculated in volt_func
    :param n:   number of samples
    :param vo:  The initial voltage supplied to the system


    :return v_o_l : returns and array of output voltage values that are a function of the time array
                    using the Runge Kutta method


    This function uses the fourth order Runge Kutta method method to simulate the discharging or charging of an
    RC circuit
    """
    h = (t_array[n-1] / (n-1))                                  # step size
    v_o_l = []                                               # list to hold output of the method
    v_o = vo                                        # initial voltage out is the voltage supplied initially
    v_o_l.append(vo)
    for j in range(0,n - 1):                                    # calculating the voltage out for time in the time array
        k1 = (1 / tow) * (vi[j] - v_o)
        k2 = (1 / tow) * ((vi[j])-(v_o + (h* k1) / 2))      # calculating the components for the Runge Kutta method
        k3 = (1 / tow) * ((vi[j])-(v_o + (h* k2) / 2))
        k4 = (1 / tow) * ((vi[j])-(v_o + h * k3))
        v_o = v_o + h * (k1 + 2*k2 + 2*k3 + k4) * (1 / 6)       # output voltage for each time element
        v_o_l.append(v_o)                                             # appending all values to the k list
    return (v_o_l)


def adams_bashforth(t_array, vi, tow , n, vo):
    """

    :param t_array:   Array of time values
    :param vi:  Array of input voltages
    :param tow: Time constant calculated in volt_func
    :param n:   number of samples
    :param vo:  The initial voltage supplied to the system


    :return x : returns and array of output voltage values that are a function of the time array using the a
    dam bashforth method

    This function uses the fourth order adam bashforth method to simulate the discharging or charging of an
    RC circuit
    """

    h = (t_array[n - 1] / (n-1))                            # step size
    v_out = vo                                          # Initial voltage out
    x = list()                                          # list to hold output of the method
    x.append(v_out)                                     # appending first voltage out to list
    x1 = v_out + (1 / tow) * (vi[1]-v_out) * h               # using Euler step to calculate second voltage out
    x.append(x1)                                        # appending second voltage out to list

    for i in range(1, 3):                              # calculate third and fourth voltage out using Euler step
        x_n1 = x[i-1] + 2*(1/tow)*(vi[i]-x[i])*h
        x.append(x_n1)                                  # Appending the individual voltage out to list

    for k in range(3,n-1):                              # Adam Bashforth method algorithm
        f_a = (1/tow) * (vi[k] - x[k])                     # calculating individual components
        f_b = (1 / tow) * (vi[k] - x[k-1])
        f_c = (1 / tow) * (vi[k] - x[k-2])
        f_d = (1 / tow) * (vi[k] - x[k-3])
        x_next = x[k] + (h / 24) * (55 * f_a - 59 * f_b + 37 * f_c - 9 * f_d) # joining all components in equation
        x.append(x_next)

    return x


def analytic_sol(t_array, v, tow):
    """
    :param t_array: time array that is used to calculate the analytic solution
    :param v: Input voltage for the circuit
    :param tow: time constant, which is resistance multiplied by capacitance


    :return v_out: returns array of voltage_out

    This function uses the analytic solution to the RC circuit and outputs an array of voltages that are a
    function of time.
    """
    v_out = []    # Creating a empty list to store the voltages

    for i in t_array:                                       # loop through the time array
        k = v * (np.exp(-(i / tow)))                        # Compute voltages using analytic solution
        v_out.append(k)                                     # Append the voltages to the list
    return v_out


def method(m , t_array, vi, tow, n, vo):
    """
    :param m: Either 'r' which stands for Runge Kutta method , or 'a' which is for the Adam Bashforth Method
    :param t_array:   Array of time values
    :param vi:  Array of input voltages
    :param tow : time constant
    :param n:   number of samples

    :return array: returns and array of values that are a function of the time array using the adam bashforth method or
    the Runge Kutta method depending on the input m, which can either be, 'r' for Runge method, or 'a' for adam method

    This Function acts as a switch, allowing for the user to decide what method: Runge Kutta or Adam Bashforth to use
    to solve the RC first oder ODE.
    """

    if m == 'r':        # if the input m is 'r' then use the Runge Kutta Method to evaluate the ODE
        array = runge_kutta(t_array, vi, tow, n, vo) # calls on pre-defined function
        return (array)
    if m == 'a':        # if the input m is 'a' then use the Adam Bashforth Method to evaluate the ODE
        array= adams_bashforth(t_array, vi, tow, n, vo) # calls on pre-defined function
        return(array)


# Define the required variables
time = 7                                            # defining the maximum value for time
voltage_in = 5                                      # Maximum input voltage
sample_num = (1) + 1000                            # Number samples, controls the step size as well,add one to fix num
res = 1000                                          # Value for the resistance  in ohms
cap = 1 * 10 ** (-3)                                # value for the capacitance in farads


# Following line calls the function to obtain the input voltage depending on what function we want,
# time array and time constant for use later on
voltage_in_array_dis, time_array_dis, tow_dis = volt_func('discharge', voltage_in, res, cap, sample_num, time)
voltage_in_array_sq1, time_array_sq1, tow_sq1 = volt_func('square1', voltage_in, res, cap, sample_num, time)
voltage_in_array_sq2, time_array_sq2, tow_sq2 = volt_func('square2', voltage_in, res, cap, sample_num, time)

# Following code plots the graphs of; Runge Kutta, Adam Bashforth and the Analytical method for no input voltage
# after t = 0
plt.figure(1)
runge_array = method('r', time_array_dis, voltage_in_array_dis, tow_dis, sample_num, voltage_in)
adam_array = method('a', time_array_dis, voltage_in_array_dis, tow_dis, sample_num , voltage_in)
analytic_array = analytic_sol(time_array_dis, voltage_in, tow_dis)


plt.xlabel("Time(s)")
plt.ylabel("V$_{out}$(v)")
plt.plot(time_array_dis, runge_array, label = 'Runge-Kutta 4$^{th}$ Order', linewidth = 9, color = 'orange')
plt.plot(time_array_dis, adam_array, label = 'Adam Bashforth 4$^{th}$ Order', linewidth = 6, color = 'm' )
plt.plot(time_array_dis, analytic_array, label = 'Analytical solution, V = V$_0$e$^{-t/RC}$', color = 'c')
plt.legend()
plt.show()

# Following code plots the graph of; Runge Kutta for the square wave input square1, which has a period of 2*r*c
runge_array_sq1 = method('r', time_array_sq1, voltage_in_array_sq1, tow_sq1, sample_num, voltage_in )

plt.xlabel("Time(s)")
plt.ylabel("V$_{out}$(v) ")
plt.plot(time_array_sq1, runge_array_sq1, label = 'Runge-Kutta V$_{out}$ 4$^{th}$ Order')
plt.plot(time_array_sq1, voltage_in_array_sq1, label = 'Driving Voltage V$_{in}$, T = 2rc' )
plt.legend(loc = 4)
plt.show()

# Following code plots the graph of; Runge Kutta for the square wave input square2, which has a period of r * c /2
runge_array_sq2 = method('r', time_array_sq2, voltage_in_array_sq2, tow_sq2, sample_num, voltage_in)

plt.xlabel("Time(s)")
plt.ylabel("V$_{out}$(v) ")
plt.plot(time_array_sq2, runge_array_sq2, label = 'Runge-Kutta V$_{out}$ 4$^{th}$ Order')
plt.plot(time_array_sq2, voltage_in_array_sq2, label = 'Driving Voltage V$_{in}$, T = rc/2')
plt.legend()
plt.show()

# Following code plots the graph of; Runge Kutta for the discharge input, for different step sizes
sample_num1 = 1 + 500
sample_num2 = 1 + 1000
sample_num3 = 1 + 2000

# Call the voltage 
voltage_in_array_1, time_array_1, tow_dis_1 = volt_func('discharge', voltage_in, res, cap, sample_num1, time)
voltage_in_array_2, time_array_2, tow_dis_2 = volt_func('discharge', voltage_in, res, cap, sample_num2, time)
voltage_in_array_3, time_array_3, tow_dis_3 = volt_func('discharge', voltage_in, res, cap, sample_num3, time)
runge_array_1 = method('r', time_array_1,voltage_in_array_1, tow_dis_1, sample_num1, voltage_in)
runge_array_2 = method('r', time_array_2, voltage_in_array_2, tow_dis_1, sample_num2, voltage_in)
runge_array_3 = method('r', time_array_3, voltage_in_array_3, tow_dis_1, sample_num3, voltage_in)

analytic_array_1 = analytic_sol(time_array_1, voltage_in, tow_dis_1)
analytic_array_2 = analytic_sol(time_array_2, voltage_in, tow_dis_2)
analytic_array_3 = analytic_sol(time_array_3, voltage_in, tow_dis_3)


# Calculating error for each line with respect to the analytical solution
difference1 = [abs(el - el2) for el, el2 in zip(runge_array_1 , analytic_array_1)]
difference2 = [abs(el - el2) for el, el2 in zip(runge_array_2 , analytic_array_2)]
difference3 = [abs(el - el2) for el, el2 in zip(runge_array_3 , analytic_array_3)]


mean1 = "%.3g" % np.mean(difference1)
mean2 = "%.3g" % np.mean(difference2)
mean3 = "%.3g" % np.mean(difference3)
print(mean1)
print(mean2)
print(mean3)


plt.xlabel("Time(s)")
plt.ylabel("V$_{out}$(v) ")

plt.plot(time_array_1, runge_array_1, label ='Runge-Kutta V$_{out}$, step size : '
                                + str(time/(sample_num1 - 1)) + ' Error : ' + str(mean1) , lw = 12, color='lightblue')
plt.plot(time_array_2, runge_array_2, label = 'Runge-Kutta V$_{out}$ , step size = '
                                + str(time/(sample_num2 - 1)) + ' Error : ' + str(mean2), lw = 9, color ='red' )
plt.plot(time_array_3, runge_array_3, label = 'Runge-Kutta V$_{out}$, step size = '
                                + str(time/(sample_num3 - 1)) + ' Error : ' + str(mean3) , lw = 6, color ='gold')
plt.plot(time_array_dis, analytic_array, label = 'Analytical solution,V = V$_0$e$^{-t/RC}$', lw = 2, color = 'b')

plt.legend()
plt.show()




