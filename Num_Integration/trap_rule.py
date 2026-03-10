def f(x): # Define the function f(x)
    pass

def simple_trap_rule(a,b,N): 
    # Simple Trapezoidal Rule on f(x) with N subintervals, 
    # approximating the integral of f from a to b
    h = (b-a)/N
    sum = 0.5*f(a)
    current = a
    for i in range(N-1):
        current = current+h
        sum += f(current)
    sum += 0.5*f(b)
    sum = sum*h
    return sum