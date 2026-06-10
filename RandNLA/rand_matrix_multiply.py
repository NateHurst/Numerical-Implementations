# Random Matrix Multiply Algorithm for approximating matrix multiplation AB.
# All matrices are numpy arrays

from math import ceil
import numpy as np

def Frob_norm_squared(A):
    A = A**2
    return np.sum(A)

def compute_probs(A):
    A_frob = Frob_norm_squared(A)
    A = A**2
    return A.sum(axis = 0)/(A_frob)

def rand_mat_multiply(A,B,epsilon): # epsilon^2 * frob(A)^2 * frob(B)^2 is the desired bound on the variance frobenius norm of AB-CR
    C_cols = []
    R_rows = []
    s = ceil(1/(epsilon**2))
    if s > len(A[0]):
        raise SystemExit("Desired epsilon bound too small, no benefit is gained from rand_matrix_multiply. (Would have to take more samples than there are columns of A or rows of B)")
    print(f"We will take {s} samples.")
    probs = compute_probs(A)
    for k in range(s):
        index = np.random.choice(len(A[0]),p = probs)
        scale = 1/(np.sqrt(s*probs[index]))
        C_cols.append(scale*A[:,index])
        R_rows.append(scale*B[index,:])
    C = np.column_stack(C_cols)
    R = np.vstack(R_rows)
    return np.dot(C,R)

# A is a n x d matix
# B is a d x m
n = 1000
d = 5000
m = 2100
A = np.random.rand(n,d)
B = np.random.rand(d,m)

epsilon = 0.01
bound = (epsilon**2)*(Frob_norm_squared(A))*(Frob_norm_squared(B))
result = np.dot(A,B)
actual = Frob_norm_squared(rand_mat_multiply(A,B,epsilon)-result)

print(f"Actual difference in Frob norm squared ||AB-CR||_F^2: {actual}")
print(f"Theoretical bound: {bound}")
print(f"Relative error ||AB-CR||_F^2 / ||AB||_F^2 : {actual/Frob_norm_squared(result)}")