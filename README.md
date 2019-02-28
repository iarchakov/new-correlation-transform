# new-correlation-transform
A New Parametrization of Correlation Matrices

Implementation of the new transformation for correlation matrices introduced in the paper "A New Parametrization of Correlation Matrices" by Archakov and Hansen (2018). The program consists of two functions. The first function, **direct_mapping_mat**, performs a direct tranform of a non-singular correlation matrix into an unconstrained real vector. The second function, **inverse_mapping_vec**, provides an inverse mapping that transforms a real vector of proper dimensionality to a unique correlation matrix that corresponds to a given vector.

Function **direct_mapping_mat** takes a correlation matrix (square, positive definite and having ones on the main diagonal) in a format of 2D numpy array as an argument. It returns a real vector of dimension *n(n-1)/2* (where *n* is a size of input matrix) in a format of 1D numpy array.

Function **inverse_mapping_vec** takes two arguments. The first argument is a real vector of size *n(n-1)/2* (for some integer *n>1*) in a format of 1D numpy array or list. The second argument is a tolerance value which has to be specified within interval *(0,1)*. The tolerance value controls the precision of convergence of the inverse mapping algorithm. The smaller value implies a more accurate convergence to the appropriate correlation matrix. The function returns the corresponding correlation matrix in a format of 2D numpy array and the number of iterations required by the algorithm to achieve a specified tolerance of convergence.
