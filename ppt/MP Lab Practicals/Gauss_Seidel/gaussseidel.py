import numpy as np
from scipy.linalg import solve
def gaussseidel(A, B):
    row, col = np.shape(A)
    if row == col:
        n = 10000
        x = B/(np.diagonal(A))
        inbuilt = solve(A,B)
        for i in range(1, n):
            x_new = np.zeros_like(x)
            print(x)
            for i in range(row):
                x_new[i] = (B[i] - np.dot(A[i, :i], x_new[:i]) - np.dot(A[i, i + 1 :], x[i + 1 :])) / A[i, i]
            if np.allclose(x, x_new, rtol=1e-8):
                break
            else:
                x = x_new
        print("Solution: ",x)
        error = np.dot(A, x) - B
        print("Inbuilt Solution",solve(A,B))
        print("Error: ",error)
    else:
        print("Matrix is not square matrix")


if __name__ == "__main__":
    A = np.array([[8, 3, -3], [-2, -8, 5], [3, 5, 10]])
    # initialize the RHS vector
    B = np.array([14,5,-8])
    # Find diagonal coefficients
    diag = np.diag(np.abs(A)) 
    # Find row sum without diagonal
    off_diag = np.sum(np.abs(A), axis=1) - diag
    if np.all(diag > off_diag):
        print('matrix is diagonally dominant')
    else:
        print('NOT diagonally dominant')
    gaussseidel(A,B)
    
   
