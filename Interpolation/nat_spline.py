def Tridagonal_solve(a,b,c,d): # Tridiagonal solver for Ax = d, where a is the diagonal of A, b is the 
    U = []                     # above-diagonal, and c is the below-diagonal. Time Complexity = O(n), n = length of diagonal.
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

def find_natspline(X,F): # X = [x_0,...,x_n], F = [f_0,...,f_n]
    a = []
    b = []
    d = []
    for j in range(1,len(X)-1):
        h1 = X[j+1]-X[j]
        h0 = X[j]-X[j-1]
        a.append(2*(h1+h0))
        if j != len(X)-2:
            b.append(h1)
        d.append((-6/h0)*(F[j]-F[j-1])+ (6/h1)*(F[j+1]-F[j]))
    M = Tridagonal_solve(a,b,b,d)
    M.insert(0,0)
    M.append(0)
    Spline = []
    for i in range(len(X)-1):
        h = X[i+1]-X[i]
        A = (1/(6*h))*(M[i+1]-M[i])
        B = 0.5*M[i]
        C = (1/h)*(F[i+1]-F[i])-(1/6)*h*(M[i+1]+2*M[i])
        D = F[i]
        Spline.append([X[i],A,B,C,D])
    return Spline

def evaluate_spline(Spline,x): # Evaluate at a point x.
    for j in range(len(Spline)):
        if x == Spline[0][0]:
            return Spline[0][1]*((x-Spline[0][0])**3)+Spline[0][2]*((x-Spline[0][0])**2)+Spline[0][3]*((x-Spline[0][0]))+Spline[0][4]
        if x <= Spline[j][0] and x >= Spline[j-1][0]:
            return Spline[j-1][1]*((x-Spline[j-1][0])**3)+Spline[j-1][2]*((x-Spline[j-1][0])**2)+Spline[j-1][3]*((x-Spline[j-1][0]))+Spline[j-1][4]
    return Spline[j][1]*((x-Spline[j][0])**3)+Spline[j][2]*((x-Spline[j][0])**2)+Spline[j][3]*((x-Spline[j][0]))+Spline[j][4]