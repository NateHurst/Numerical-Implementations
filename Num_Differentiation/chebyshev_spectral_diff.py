# Problem 7.7.7
# Euler-Maclaurin modified quadrature

from math import exp,cos,pi,sqrt,sin
import numpy as np

def f(x):
    return exp(x)

def simple_trap_rule(a,b,N):
    h = (b-a)/N
    sum = 0.5*f(a)
    current = a
    for i in range(N-1):
        current = current+h
        sum += f(current)
    sum += 0.5*f(b)
    sum = sum*h
    return sum

def discrete_cos_transform(f_vals):
    output = np.zeros(len(f_vals),dtype=float)
    for k in range(len(f_vals)):
        sum = 0.5*(f_vals[0]+f_vals[-1]*cos(k*pi))
        for j in range(1,len(f_vals)-1):
            sum += f_vals[j]*cos((j*k*pi)/(len(f_vals)-1))
        output[k] = sum
    return (2/(len(f_vals)-1))*output

def fast_cheby_diff(f,a,b,N,func_vals = None): # Fast Chebyshev differentiation in [-1,1] maped to [a,b]
    nodes = [cos((j*pi)/N) for j in range(N+1)]
    if func_vals is None:
        func_vals = [f((a+b)/2 + node*((b-a)/2)) for node in nodes]
    func_vals = discrete_cos_transform(func_vals)
    first = ((N**2)*func_vals[len(func_vals)-1])/2
    for k in range(len(func_vals)-1):
        first += (k**2)*func_vals[k]
    last = 0.5*((-1)**(N+1))*N*N*func_vals[-1]
    for k in range(len(func_vals)-1):
        last += ((-1)**(k+1))*(k**2)*func_vals[k]
    modified = [0]
    for k in range(1,len(func_vals)-1):
        modified.append((-1)*k*func_vals[k])
    modified.append(0)
    pi_prime_vals = np.zeros(len(func_vals),dtype=float)
    for r in range(1,len(func_vals)-1):
        for p in range(len(modified)):
            pi_prime_vals[r] += modified[p]*(sin(p*r*pi/N))
    for i in range(1,len(func_vals)-1):
        func_vals[i] = (pi_prime_vals[i])/((-1)*(sqrt(1-nodes[i]**2)))
    func_vals[0]= first
    func_vals[N]= last
    return (func_vals*(2/(b-a)))


def derivs(f,num_derivs,a,b,N):
    output = fast_cheby_diff(f,a,b,N)
    for k in range(1,num_derivs):
        output = fast_cheby_diff(f,a,b,N,output)
    return output
    
def Quad(f,a,b,N):
    h = (b-a)/N
    out = simple_trap_rule(a,b,N)
    der1 = derivs(f,1,a,b,8)
    der3 = derivs(f,3,a,b,8)
    out += ((-1)/12)*(h**2)*(der1[0]-der1[-1])
    out += (1/(24*30))*(h**4)*(der3[0]-der3[-1])
    return out

def Composite_mod_trap_rule(f,a,b,N):
    h = (b-a)/N
    sum = 0.5*(f(a)+f(b))
    for k in range(1,N):
        sum += f(a+k*h)
    sum *= h
    der = derivs(f,1,a,b,8)
    sum -= ((h**2)/12)*(der[0]-der[-1])
    return sum

a = 0
b = 1
approx_quad = [Quad(f,a,b,n) for n in [2,4,8,16]]
error_quad = [abs(app-(exp(1)-exp(0))) for app in approx_quad]
approx_mod_trap = [Composite_mod_trap_rule(f,a,b,n) for n in [2,4,8,16,32,64]]
error_mod_trap = [abs(app-(exp(1)-exp(0))) for app in approx_mod_trap]

for i in range(1,len(error_quad)):
    print(f"Approximate convergence rate (Quadrature) Err_N/Err_(N/2)) is {error_quad[i-1]/error_quad[i]}")

print("For quad:")
print(approx_quad)
print(error_quad)
print("For mod trap:")
print(approx_mod_trap)
print(error_mod_trap)