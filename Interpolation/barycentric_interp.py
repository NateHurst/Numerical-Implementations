# Input: Evaluation point x, barycentric weights lambda_j, data (x_j,f_j)
# Output: Interpolating polynomial using barycentric formula of interpolation
# Time Complexity: O(n^2) one time cost to find the weights (precomputed), 
# O(n) to evaluate the interpolating polynomial at a point x.

from math import cos, exp

def f(x): # Define the function f(x) if using it to find data points
    pass

def compute_weights(Data): # Data is an array of arrays [[x_0,f_0],[x_1, f_1],...,[x_n,f_n]]
    weights = {}
    for j in range(len(Data)):
        omega_prime = 1
        for k in range(len(Data)):
            if k !=j:
                omega_prime *= Data[j][0] - Data[k][0]
        key = j+1
        weights[key] = 1/(omega_prime)
    return weights

def interpolate(x,Data,weights):
    if [x,f(x)] in Data:
        return f(x)
    else:
        sum1 =0
        sum2 =0
        for i in range(len(Data)):
            sum1 += (weights[i+1]/(x-Data[i][0]))*Data[i][1]
            sum2 += weights[i+1]/(x-Data[i][0])
        return (sum1/sum2)