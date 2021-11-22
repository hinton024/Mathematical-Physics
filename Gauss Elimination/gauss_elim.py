import numpy as np     
def elimination(a,b):
    row, col = np.shape(a)
    if row == col:        #condition to check square matrix
        n = len(b)
        x  = np.zeros(n, float)
        for p in range(n-1):       #index the pivots
            if a[p,p] < 1e-12: #putting threshold condition to check pivot not to be zero
                for i in range(p+1,n):  #index elements under pivot
                    if np.fabs(a[i,p]) > np.fabs(a[p,p]):  #non zero rows under pivot element
                        a[[p,i]] = a[[i,p]]   #interchange rows and columns 
                        b[[p,i]] = b[[i,p]]
                        break
            for i in range(p+1,n):      
                if a[i,p] == 0:    #skip the row 
                    continue
                factor = a[p,p]/a[i,p]
                for j in range(p,n):
                    a[i,j] = a[p,j] - a[i,j]*factor
                b[i] = b[p] - b[i]*factor     #elimination
        '''check rank of matrix to determine consistency'''
        rank_a = np.linalg.matrix_rank(a)   
        aug_matrix = np.column_stack((a,b.T))
        aug_matrix_rank = np.linalg.matrix_rank(aug_matrix)
        if rank_a == aug_matrix_rank:
            if aug_matrix_rank == np.shape(b)[0]:
                print("System of equations has unique solution")
                print("Augmented Matrix \n" ,np.column_stack((a,b.T)))
                '''Back Substitution'''
                x[n-1] = b[n-1] / a[n-1,n-1]
                for i in range(n-2,-1, -1):
                    sum_ax = 0
                    for j in range(i+1, n):
                        sum_ax += a[i, j]* x[j]
                    x[i] = (b[i] - sum_ax) / a[i,i]

                print('The solution of the system: ')
                print(x)
                numpy_solution = np.linalg.solve(a, b)
                print("Solution by inbuilt function", numpy_solution)
            elif aug_matrix_rank < np.shape(b)[0]:
                print("system has infinitely many solutions")
        elif rank_a < aug_matrix_rank:
            print("System is inconsistent")
    else:
        print("Number of rows should be equal to number of columns")


if __name__ == "__main__":
    print("E.g. 1")
    a = np.array([[1,-2,1],
                [2,-5,4],
                [1,-4,6]], float)
    b = np.array([5,-3,10], float)
    elimination(a,b)

    print("E.g. 2")
    a = np.array([[1,-2,1,-1,1],
                [2,-5,4,1,-1],
                [1,-4,6,2,-1]], float)
    b = np.array([5,-3,10], float)
    elimination(a,b)
    print("E.g. 3")
    a = np.array([[1,-5,4],
                [1,-5,3],
                [2,-10,13]], float)
    b = np.array([3,6,5], float)
    elimination(a,b)

    print("E.g. 4")
    a = np.array([[12,-3],
                [-12,3]], float)
    b = np.array([6,-6], float)
    elimination(a,b)
    
    