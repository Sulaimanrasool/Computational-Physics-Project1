
# This File answers question 1 of the Assignment

def near(x):
    """
    :param x: A floating number

    :return up_val: Returns nearest representable real number higher than x
    :return low_val: Returns nearest representable real number lower than x
    :return ud: Returns upper difference within which numbers are rounded to x
    :return ld : Returns the lower difference within which numbers are rounded to x
    :return frac: Returns fractional rounding range

    This function takes an input value x and returns the parameters above
    """
    upper = (x + x/2)                       # Initial upper value, we iterate down from here
    lower = (x - x/2)                       # Initial lower value, we iterate up from here
    u_m = list()                            # list to store upper value
    u_l = list()                            # # list to store lower value
    counter1 = 0                            # required for while loop check
    counter2 = 0
    while upper != x:                       # Iteration to get to the upper nearest representable number
        counter1 += 1                       # increment counter for break loop check
        upper = 1/2 * (upper + x)           # Iterate the value closer to x
        if upper != x:
            u_m.append(upper)
        if counter1 > 2:                    # If statement prevents loop extending forever
            if u_m[-1] == u_m[-2]:
                break

    while lower != x:                     # Iteration to get to the lower nearest representable number
        counter2 += 1
        lower = 1 / 2 * (lower + x)         # Iterate closer to x
        if lower != x:
            u_l.append(lower)

        if counter2 > 2:                    # If statement prevents loop extending forever
            if u_l[-1] == u_l[-2]:
                break

    up_val = u_m[-1]                        # Obtains last value in the list for upper
    low_val = u_l[-1]                       # Obtains last value in the list for lower
    ud = up_val - x                         # Get the upper difference
    ld = x - low_val                        # Get the lower difference
    frac = (((ud/2) + (ld/2))/x)

    return up_val, low_val, ud, ld, frac


# The following code returns outputs for the input value x = 0.25 and prints the fractional rounding range of 0.25
upper_value, lower_value, upper_diff, lower_diff, frac_range = near(0.25)
print("The Upper value of 0.25 is " + str(upper_value))
print("The lower value of 0.25 is " + str(lower_value))
print("The Fractional rounding range of 0.25 is: " + "" + str(frac_range))

print(" ")

# Following code uses the two closest values to 0.25, runs the function near, obtaining closest values to them in turn
# and returns the fractional rounding range of the  nearest representable real number higher than 0.25 and returns the
# fractional rounding range of the  nearest representable real number lower than 0.25
upper_value_1, lower_value_1, upper_diff_1, lower_diff_1, frac_range_1 = near(upper_value)
upper_value_2, lower_value_2, upper_diff_2, lower_diff_2, frac_range_2 = near(lower_value)


print("The Upper value of " + str(upper_value) + " is: " + str(upper_value_1))
print("The lower value of " + str(upper_value) + " is: " + str(lower_value_1))
print("The Fractional rounding range of " + str(upper_value) + " is: " + "" + str(frac_range_1))
print("")


print("The Upper value of " + str(lower_value) + " is: " + str(upper_value_2))
print("The lower value of " + str(lower_value) + " is: " + str(lower_value_2))
print("The Fractional rounding range of " + str(lower_value) + " is: " + "" + str(frac_range_2))

# Fractional range change because you have to change the mantissa and the exponent

