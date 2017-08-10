int dbscan(long **all_points,
           long **dist_matrix,
           long *cluster_assignments, 
           int num_points,
           int num_dimensions,
           int alphabet_size,
           int epsilon_cutoff,
           int min_neighbors,
           int dist_function);
