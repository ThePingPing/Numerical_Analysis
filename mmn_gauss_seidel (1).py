def gauss_seidel(A, b, x, tol, max_iter):
    k = 1
    n = len(b)
    xo = x[:]
    while k <= max_iter:
        for i in range(n):
            if A[i][i] == 0:
                print("matrix not suitable for this convergence method")
                return
            x[i] = (1/A[i][i]) * \
                   (b[i] - sum([A[i][j] * x[j] for j in range(i)]) - sum([A[i][j] * xo[j] for j in range(i+1, n)]))
        if max([abs(x[i] - xo[i]) for i in range(n)]) < tol:
            print(f"number of iterations: {k}")
            return x
        k += 1
        xo = x[:]


if __name__ == "__main__":
    A = [[2,-1,1],[2,2,2],[-1,-1,2]]
    b = [-1,4,-5]
    x = [0,0,0]
    tol = 0.00001
    max_iter = 1000
    print(gauss_seidel(A, b, x, tol, max_iter))



