cimport numpy as np
np.import_array()

from libc.stdlib cimport malloc, free

# cdefine the signature of our c function
cdef extern from "dbscan.h":
    int dbscan(long **all_points,
               long **dist_matrix,
               long *cluster_assignments,
               int num_points,
               int num_dimensions,
               int alphabet_size,
               int epsilon_cutoff,
               int min_neighbors,
               int dist_function)

# define our simple wrap function
def run_dbscan(np.ndarray[long, ndim=2, mode="c"] all_points not None, 
                np.ndarray[long, ndim=2, mode="c"] dist_matrix not None,
                np.ndarray[long, ndim=1, mode="c"] cluster_assignments not None,
                num_points,
                num_dimensions,
                alphabet_size,
                epsilon, 
                min_neighbors,
                dist_function):

    cdef:
        long **_all_points
        long **_dist_matrix

    _all_points = <long **> malloc(sizeof(long *) * num_points)
    _dist_matrix = <long **> malloc(sizeof(long *) * alphabet_size)
    
    for i in range(num_points):
        _all_points[i] = <long *> np.PyArray_DATA(all_points[i,:])

    for i in range(alphabet_size):
        _dist_matrix[i] = <long *> np.PyArray_DATA(dist_matrix[i,:])

    status = dbscan(_all_points,
                    _dist_matrix,
                    <long*> np.PyArray_DATA(cluster_assignments),
                    <int> num_points,
                    <int> num_dimensions,
                    <int> alphabet_size,
                    <int> epsilon,
                    <int> min_neighbors,
                    <int> dist_function)

    free(_all_points)
    free(_dist_matrix)

    return status
