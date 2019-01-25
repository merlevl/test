import numpy as np


####################################################################################################
# Exercise 1: Interpolation

def lagrange_interpolation(x: np.ndarray, y: np.ndarray) -> (np.poly1d, list):
    """
    Generate Lagrange interpolation polynomial.

    Arguments:
    x: x-values of interpolation points
    y: y-values of interpolation points

    Return:
    polynomial: polynomial as np.poly1d object
    base_functions: list of base polynomials
    """

    assert (x.size == y.size)

    polynomial = np.poly1d(0)
    lagrange = np.poly1d(1)
    base_functions = []
    # TODO: Generate Lagrange base polynomials and interpolation polynomial
    for i in range(x.size):
        for j in range(x.size):
            rechen = np.poly1d([1, -1 * x[j] ])
            if (i < j) :
                lagrange = lagrange * rechen / (x[i] - x[j])
            elif (i > j):
                lagrange = lagrange * rechen / (x[i] - x[j])
        base_functions.append(lagrange)
        polynomial += (lagrange * y[i])
        lagrange = np.poly1d(1)
    return polynomial, base_functions

def hermite_cubic_interpolation(x: np.ndarray, y: np.ndarray, yp: np.ndarray) -> list:
    """
    Compute hermite cubic interpolation spline

    Arguments:
    x: x-values of interpolation points
    y: y-values of interpolation points
    yp: derivative values of interpolation points

    Returns:
    spline: list of np.poly1d objects, each interpolating the function between two adjacent points
    """

    assert (x.size == y.size == yp.size)

    spline = []
    # TODO compute piecewise interpolating cubic polynomials
    M = np.array([[1.0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 1, 0, 0]])
    F = np.empty(4)
    for i in range(x.size - 1):
        M[0][1] = x[i]
        M[0][2] = x[i] ** 2
        M[0][3] = x[i] ** 3
        M[1][1] = x[i + 1]
        M[1][2] = x[i + 1] ** 2
        M[1][3] = x[i + 1] ** 3
        M[2][2] = x[i] *  2
        M[2][3] = (x[i] ** 2) *3
        M[3][2] = x[i + 1] *  2
        M[3][3] = (x[i + 1] ** 2) *3
        
        F[0]    = y[i]
        F[1]    = y[i + 1]
        F[2]    = yp[i]
        F[3]    = yp[i + 1]
        
        spline.append(np.poly1d(np.linalg.solve(M, F)[::-1]))
    return spline



####################################################################################################
# Exercise 2: Animation

def natural_cubic_interpolation(x: np.ndarray, y: np.ndarray) -> list:
    """
    Intepolate the given function using a spline with natural boundary conditions.

    Arguments:
    x: x-values of interpolation points
    y: y-values of interpolation points

    Return:
    spline: list of np.poly1d objects, each interpolating the function between two adjacent points
    """

    assert (x.size == y.size)
    # TODO construct linear system with natural boundary conditions
    M = np.zeros((4 * (x.size - 1), 4 * (x.size -1 )), dtype = np.float64)
    F = np.zeros(4 * (x.size - 1))
    for i in range(x.size - 1): 
        M[4 * i][4 * i] = 1
        M[4 * i][4 * i + 1] = x[i]
        M[4 * i][4 * i + 2] = x[i] ** 2
        M[4 * i][4 * i + 3] = x[i] ** 3
        
        M[4 * i + 1][4 * i] = 1
        M[4 * i + 1][4 * i + 1] = x[i + 1]
        M[4 * i + 1][4 * i + 2] = x[i + 1] ** 2
        M[4 * i + 1][4 * i + 3] = x[i + 1] ** 3
        
        F[i * 4] = y[i]
        F[i * 4 + 1] = y[i + 1]
        
    for i in range(x.size - 2):    
        
        M[4 * i + 2][4 * i + 1] = 1
        M[4 * i + 2][4 * i + 2] = x[i + 1] * 2
        M[4 * i + 2][4 * i + 3] = x[i + 1] ** 2 *3
        
        M[4 * i + 2][4 * i + 5] = -1
        M[4 * i + 2][4 * i + 6] = x[i + 1] * -2
        M[4 * i + 2][4 * i + 7] = x[i + 1] ** 2 * -3
        
        M[4 * i + 3][4 * i + 2] = 2
        M[4 * i + 3][4 * i + 3] = x[i + 1] * 6
        
        M[4 * i + 3][4 * i + 6] = -2
        M[4 * i + 3][4 * i + 7] = x[i + 1] * -6
        
    M[4 * x.size - 6][2] = 2
    M[4 * x.size - 6][3] = 6 * x[0]
    
    M[4 * x.size - 5][4 * x.size - 6] = 2
    M[4 * x.size - 5][4 * x.size - 5] = 6 * x[x.size - 1]
    
    # TODO solve linear system for the coefficients of the spline
    print("\n", x.size, "\n")
    for i in range(x.size):
        print(x[i], "\n")
    print(M)
    c = np.linalg.solve(M, F)

    spline = []
    # TODO extract local interpolation coefficients from solution
    for i in range(x.size - 1):
        indices = [4 * i + 3, 4 * i + 2, 4 * i + 1, 4 * i]
        spline.append(np.poly1d(np.take(c, indices)))

    return spline


def periodic_cubic_interpolation(x: np.ndarray, y: np.ndarray) -> list:
    """
    Interpolate the given function with a cubic spline and periodic boundary conditions.

    Arguments:
    x: x-values of interpolation points
    y: y-values of interpolation points

    Return:
    spline: list of np.poly1d objects, each interpolating the function between two adjacent points
    """

    assert (x.size == y.size)
    # TODO: construct linear system with periodic boundary conditions
    M = np.zeros((4 * (x.size - 1), 4 * (x.size -1 )), dtype = np.float64)
    F = np.zeros(4 * (x.size - 1))
    for i in range(x.size - 1): 
        M[4 * i][4 * i] = 1
        M[4 * i][4 * i + 1] = x[i]
        M[4 * i][4 * i + 2] = x[i] ** 2
        M[4 * i][4 * i + 3] = x[i] ** 3
        
        M[4 * i + 1][4 * i] = 1
        M[4 * i + 1][4 * i + 1] = x[i + 1]
        M[4 * i + 1][4 * i + 2] = x[i + 1] ** 2
        M[4 * i + 1][4 * i + 3] = x[i + 1] ** 3
        
        F[i * 4] = y[i]
        F[i * 4 + 1] = y[i + 1]
        
    for i in range(x.size - 2):    
        
        M[4 * i + 2][4 * i + 1] = 1
        M[4 * i + 2][4 * i + 2] = x[i + 1] * 2
        M[4 * i + 2][4 * i + 3] = x[i + 1] ** 2 *3
        
        M[4 * i + 2][4 * i + 5] = -1
        M[4 * i + 2][4 * i + 6] = x[i + 1] * -2
        M[4 * i + 2][4 * i + 7] = x[i + 1] ** 2 * -3
        
        M[4 * i + 3][4 * i + 2] = 2
        M[4 * i + 3][4 * i + 3] = x[i + 1] * 6
        
        M[4 * i + 3][4 * i + 6] = -2
        M[4 * i + 3][4 * i + 7] = x[i + 1] * -6
    
    M[4 * x.size - 6][1] = 1
    M[4 * x.size - 6][2] = 2 * x[0]
    M[4 * x.size - 6][3] = 3 * (x[0] ** 2)
    
    M[4 * x.size - 6][4 * x.size - 7] = -1
    M[4 * x.size - 6][4 * x.size - 6] = -2 * x[x.size - 1]
    M[4 * x.size - 6][4 * x.size - 5] = -3 * (x[x.size - 1] ** 2)
    

    M[4 * x.size - 5][2] = 2 
    M[4 * x.size - 5][3] = 6 * x[0]
    

    M[4 * x.size - 5][4 * x.size - 6] = -2
    M[4 * x.size - 5][4 * x.size - 5] = -6 * x[x.size - 1]
    # TODO solve linear system for the coefficients of the spline
    print("\n", x.size, "\n")
    for i in range(x.size):
        print(x[i], "\n")
    print(M)
    c = np.linalg.solve(M, F)
    spline = []

    # TODO extract local interpolation coefficients from solution
    for i in range(x.size - 1):
        indices = [4 * i + 3, 4 * i + 2, 4 * i + 1, 4 * i]
        spline.append(np.poly1d(np.take(c, indices)))

    return spline


if __name__ == '__main__':

    print("All requested functions for the assignment have to be implemented in this file and uploaded to the "
          "server for the grading.\nTo test your implemented functions you can "
          "implement/run tests in the file tests.py (> python3 -v test.py [Tests.<test_function>]).")
