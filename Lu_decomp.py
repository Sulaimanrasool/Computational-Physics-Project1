import numpy as np


# This python File answers question 2 of the Assignment

def decompose(A):
    """
    This function decompose takes in a matrix A of NxN dimensions and decomposes it through clouts algorithm
    :param A: Matrix NxN
    :return L : lower diagonal matrix
    :return U : upper diagonal matrix
    :return det : determinant for the matrix A

    """

    X_dim, Y_dim = np.shape(A)                 # obtaining the dimensions of the matrix A
    L = np.zeros((X_dim, Y_dim))               # Create a matrix of size A
    U = np.zeros((X_dim, Y_dim))               # Create a matrix of size U
    det = 1                                    # Initial determinant
    for i in range(X_dim):                      # Iterate over X dimension
        for j in range(Y_dim):                  # Iterate over Y dimension
            if i <= j :                         # This loop creates the upper half of matrix U and L
                if i == j :                     # diagonal ones for matrix L
                    L[i][j] = 1
                else:
                    L[i][j] = 0
                tots = 0                         # Calculating the values of L and U in upper matrix
                for k in range(i+1):
                    tots = tots + L[i][k] * U[k][j]
                U[i][j] = A[i][j] - tots

            if i >= j :                         # This loop creates the lower half of matrix U and L
                if i ==j:                       # Creating the diagonal ones for matrix L
                    L[i][j] = 1

                else:                           # Calculating the values of L and U in lower Matrix
                    U[i][j] = 0
                    n_sum = 0
                    for r in range(j + 1):
                        n_sum = n_sum + L[i][r] * U[r][j]
                    L[i][j] = (1 / U[j][j]) * (A[i][j] - n_sum)
    for i in range(X_dim):                      # Calculates the determinant for the matrix A
        for j in range(Y_dim):
            if i == j :
                det = det * U[i][j]
    return L, U, det


def LUdecompose(A,b):
    """"
    Function requires an input Matrix A, NxN dimensions and matrix b, Nx1 dimensions
    calculates the matrix X for , A*X = b
    outputs X.
    Requires the use of function decompose
    :param A : Matrix NxN
    :param b : matrix Nx1
    :return x : matrix for the equation A*x = B

    """
    L, U, D = decompose(A)                      # Call the function decompose to breakdown matrix A
    X_dim, Y_dim = np.shape((b))                # setting required parameters
    y = np.zeros((X_dim, Y_dim))                # Creating a matrix of zeros
    x = np.zeros((X_dim, Y_dim))
    for i in range(X_dim):                      # Forward substitution
        if i == 0:
            y[i] = b[0] / L[0][0]               # Using the equation for forward sub
        else:
            tots = 0
            for k in range(i + 1):
                tots = tots + L[i][k] * y[k]    # updating tots value
            y[i] = (1 / L[i][i]) * (b[i] - tots)

    for i in range(X_dim - 1, 0 - 1, -1):       # Backward Substitution
        if i == (X_dim - 1):
            x[i] = y[i] / (U[i][i])
        else:
            tots = 0
            for k in range(X_dim - 1):
                tots = tots + U[i][k + 1] * x[k + 1]
            x[i] = (1 / U[i][i]) * (y[i] - tots)        # update the values in matrix x
    return (x)


def inverse_matrix(A,b):
    """
    This function outputs the inverse of provided matrix A.
    Requires the function LUdecompose

    :param A: matrix A, NxN
    :param b: matrix b, Nx1

    :return A_inv: Inverse of matrix A
    """
    X_dim, Y_dim = np.shape(A)                   # obtaining the dimensions of matrix A
    id = np.identity((X_dim))                    # Creating an identity matrix
    A_inv = np.zeros((X_dim , Y_dim))            # creating a place holder matrix for the inverse of matrix A
    for i in range(X_dim):                       # Calculating the column vectors for the inverse matrix
        k = id[:,i]
        b = k.reshape(X_dim, 1)
        x = LUdecompose(A,b)                      # Call my LUdecompose function to get the matrix x
        A_inv[:, i] = x.reshape(X_dim)            # Reshape my x array and replace the ith column in A_inv with x
    return A_inv


def main():
    """"
    This function prevents the print statements from being called in interpolation
    """
    # The following code is setting up the variables that will be used for the functions
    A = np.array([[3,1,0,0,0], [3,9,4,0,0], [0,8,20,10,0], [0,0,-22,31,-25], [0,0,0,-35,61]])   # Matrix A
    b = np.array([[2], [5], [-4], [8], [9]])                                                    # Eigenvalues b

    # The following code expresses the matrix A as a upper and lower diagonal matrix and computes the determinant of A
    lower_matrix, upper_matrix, determinant = decompose(A)
    print("The lower diagonal Matrix is :" )
    print(lower_matrix)
    print("")
    print("The Upper diagonal Matrix is :" )
    print(upper_matrix)
    print("")
    print("The determinant for Matrix A is : " + "" + str(determinant))
    print("")


    # The following code uses the above functions to return the matrix x for the equation A*x = b

    x_matrix = LUdecompose(A,b)
    print("The matrix x is")
    print(x_matrix)
    print("")

    # The following code will run the function inverse to fine the inverse of Matrix A
    A_inverse_matrix = inverse_matrix(A,b)
    print("The inverse of matrix A is :")
    print(A_inverse_matrix)
    print("")

    # The following code checks if the inverse of the matrix is correct, by dotting the inverse by the original matrix
    # which results in an identity matrix
    print("The dot product of the inverse matrix with the original matrix")
    I = np.dot(A_inverse_matrix, A)
    print(I)


if __name__ == "__main__":
    main()


