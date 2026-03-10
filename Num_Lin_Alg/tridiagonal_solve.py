def Tridagonal_solve(a,b,c,d): # Tridiagonal solver for Ax = d, where a is the diagonal of A, b is the 
    U = []                     # above-diagonal, and c is the below-diagonal. Time Complexity = O(n), n= length of diagonal
    L = []                     
    Y = []
    X = []
    U.append(a[0])
    for j in range(len(a)-1):
        L.append(c[j]/U[j])
        U.append(a[j+1]-L[j]*b[j])
    Y.append(d[0])
    for j in range(1,len(a)):
        Y.append(d[j]-L[j-1]*Y[j-1])
    X.append(Y[len(a)-1]/U[len(a)-1])
    for j in range(len(a)-2,-1,-1):
        X.insert(0,(Y[j]-b[j]*X[0])/U[j])
    return X