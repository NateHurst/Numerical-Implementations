def f(x): # Function whose derivative is to be approximated
    return exp(x)

def g(x,theta,num_derivs):
    return f(x+exp(1j*theta))*exp(-1j*num_derivs*theta)

def cauchy_trap_rule(g,x,a,b,N,num_derivs):
    h = (b-a)/N
    sum = 0.5*g(x,a,num_derivs)
    current = a
    for i in range(N-1):
        current = current+h
        sum += g(x,current,num_derivs)
    sum += 0.5*g(x,b,num_derivs)
    sum = sum*h
    return sum

def cauchy_diff(x,num_derivs,N): # Using Cauchy's integral forumla and trapezoidal rule to
                                 # approximate higher order derivatives of f (f must be complex
                                 # analytic, or be able to be extended as such)
    return np.real(cauchy_trap_rule(g,x,0,2*pi,N,num_derivs)*(factorial(num_derivs)/(2*pi)))