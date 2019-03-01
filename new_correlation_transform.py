# Programming implementation of the new method of unconstrained transformation
# for correlation matrices suggested in Archakov and Hansen (2018)
# ----------------------------------------------------------------------------


import numpy as np
from scipy.linalg import logm, expm, norm



### New tranformation of correlation matrices (from a non-singular correlation matrix "C"
### to a unique corresponding real vector "gamma")
def direct_mapping_mat(C):
    gamma = []

    try:
        # Check if input is of proper format: C is 2D np.array of suitable dimensions
        # and is positive-definite correlation matrix
        if not isinstance(C, np.ndarray):
            raise ValueError
        if C.ndim != 2:
            raise ValueError
        if C.shape[0] != C.shape[1]:
            raise ValueError
        if not all([np.array_equal(np.diag(C),np.ones(C.shape[0])), np.all(np.linalg.eigvals(C) > 0),
                    np.array_equal(C,C.T)]):
            raise ValueError

        # Apply matrix log-transformation to C and get off-diagonal elements
        A = logm(C)
        gamma = A[np.triu_indices(C.shape[0], 1)]

    except ValueError:
        print('Error: input is of a wrong format')

    return gamma



### Inverse tranformation (from a real vector "gamma" of proper dimensionality
### to a unique corresponding correlation matrix "C")
def inverse_mapping_vec(gamma, tol_value=1e-8):
    C = []
    iter_number = -1

    try:
        # Check if input is of proper format: gamma is of suitable length
        if not isinstance(gamma, (np.ndarray, list)):
            raise ValueError
        if isinstance(gamma, np.ndarray):
            if gamma.ndim != 1:
                raise ValueError
        n = 0.5 * (1 + np.sqrt(1 + 8 * len(gamma)))
        if not all([n.is_integer(), n > 1]):
            raise ValueError
            
        # Check if tolerance value belongs to a proper interval
        # and change it to the default value otherwise
        if not (0 < tol_value <= 1e-4):
            tol_value = 1e-8
            print('Warning: tolerance value has been changed to default')
            
        # Place elements from gamma into off-diagonal parts
        # and put zeros on the main diagonal of [n x n] symmetric matrix A
        n = int(n)
        A = np.zeros(shape=(n,n))
        A[np.triu_indices(n,1)] = gamma
        A = A + A.T

        # Read properties of the input matrix
        diag_vec = np.diag(A)                # get current diagonal
        diag_ind = np.diag_indices_from(A)   # get row and column indices of diagonal elements

        # Iterative algorithm to get the proper diagonal vector
        dist = np.sqrt(n)
        while dist > np.sqrt(n) * tol_value:
            diag_delta = np.log(np.diag(expm(A)))
            diag_vec = diag_vec - diag_delta
            A[diag_ind] = diag_vec
            dist = norm(diag_delta)
            iter_number += 1

        # Get a unique reciprocal correlation matrix
        C = expm(A)
        np.fill_diagonal(C, 1)

    except ValueError:
        print('Error: input is of a wrong format')

    return C, iter_number
