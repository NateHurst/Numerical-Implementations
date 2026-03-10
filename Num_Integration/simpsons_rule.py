def f(x): # Define the function f(x).
    pass

def S(a,b,N): # Simpsons rule to approximate the integral of f from a to b with N subintervals
              # Same computational cost as the Trapezoid rule with 2*N intervals.
    h = (b-a)/(float(2*N))
    sum1 = (f(a)+f(b))/(2.0)
    sum2 = sum1
    point = a
    for i in range(0,(2*N)-1):
        point += h
        current = f(point)
        sum1 += current
        if i % 2 == 1:
            sum2 += current
    sum1 = h*sum1
    sum2 = (2*h)*sum2
    return (4*sum1-sum2)/(3.0)