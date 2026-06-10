from math import sqrt
import numpy as np
from numpy.linalg import eig

def Euclidean_norm(v):
    sum = 0
    for i in range(len(v)):
        sum+= v[i]**2
    return sqrt(sum)

# Find absolute value of first k eigenvalues and eigenvectors of a symmetric matrix A, with N iterations

def Power_Method(A,k,N):
    dim = len(A[0])
    vectors = np.zeros((k,dim))
    for i in range(k):
        vectors[i] = np.random.rand(dim)
        vectors[i] /= Euclidean_norm(vectors[i])
    for j in range(N): 
        for i in range(k): # Apply A
            vectors[i] = np.dot(A,vectors[i])
        for r in range(k): # Re-orthogonalize and normalize
            for p in range(r):
                vectors[r] -= np.dot(vectors[r],vectors[p])/(Euclidean_norm(vectors[p])**2)*vectors[p]
            vectors[r] = vectors[r]/Euclidean_norm(vectors[r])
    eigenvalues = np.zeros(k)
    for i in range(k): # Eigenvalue approx using Rayleigh quotient
        eigenvalues[i] = (np.dot(vectors[i],np.dot(A,vectors[i])))/(np.dot(vectors[i],vectors[i]))
    return eigenvalues, vectors

N = 20
A = np.random.rand(N,N)
A = (A+A.T)/2

eigenvalues_real, eigenvectors_real = eig(A)

k= N
approx_eigenvalues, approx_eigenvectors = Power_Method(A,k,200)
approx_eigenvalues = np.sort(approx_eigenvalues)[::-1]
eigenvalues_real = np.sort(eigenvalues_real)[::-1]

print(max(abs(approx_eigenvalues-eigenvalues_real[0:k])))
print("APPROXIMATIONS ___________________________________")
print(approx_eigenvalues)
print("REAL ___________________________________")
print(eigenvalues_real[0:k])





    
