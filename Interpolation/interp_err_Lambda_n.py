# an upper bound for the interpolation error from the interpolating polynomial p_n to f given in terms of 
# The best uniform approximation to f p^* is given by ||f-p_n||_oo <= max(L(x)) ||f-p^*||_oo where L is defined below
# as the sum of the absolute value of the cardinal functions for the nodes x_0,...,x_n

def L(x,A): # A is an array of nodes [x_0,x_1,...,x_n]
    def l(x,j): # l_j (x) function
        up_prod = 1
        low_prod = 1
        for i in range(len(A)):
            if i != j:
                up_prod = up_prod*(x-A[i])
                low_prod = low_prod*(A[j]-A[i])
        return (up_prod/low_prod)
    sum = 0
    for i in range(len(A)):
        sum += abs(l(x,i))
    return sum
