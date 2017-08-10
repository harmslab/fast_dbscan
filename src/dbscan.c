#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    
    int num_points, num_dimensions, min_neighbors, alphabet_size;
    int epsilon_cutoff;

    long **all_points;
    long *cluster_assignments;

    long **dist_matrix;

    int (*dist_function)(long *, long *, void *);
    int **dl_d, *dl_da;

} config;

/*
void loadData(char file_name[],config *c) {

    int i, j;

    c->all_points = (long**)calloc(c->num_points, sizeof(long*));
    c->cluster_assignments = (long*)calloc(c->num_points, sizeof(long));
    
    for(i = 0; i < c->num_points; i++) {
        c->all_points[i] = (long*)calloc(c->num_dimensions, sizeof(long));
    }

    // -1 indicates that sequence has not yet been assigned to cluster. 
    for (i = 0; i < c->num_points; i++){
        c->cluster_assignments[i] = -1;
    }

    c->dist_matrix = malloc(c->alphabet_size*sizeof(long *)); 
    for (i = 0; i < c->alphabet_size; i++){

        c->dist_matrix[i] = malloc(c->alphabet_size*sizeof(long));
        for (j = 0; j < c->alphabet_size; j++){
            
            if (i == j) { 
                c->dist_matrix[i][j] = 0.0;
            } else { 
                c->dist_matrix[i][j] = 1.0;
            }

        }

    }
 
    FILE *fp;
    fp = fopen(file_name, "r");
    
    for (i = 0; i < c->num_points; i++) {
        for (j = 0; j < c->num_dimensions; j++) {
            fscanf(fp, "%d", &(c->all_points[i][j]));
        }
    }
    
    fclose(fp);
} */


int simple_distance(long *v1, long *v2, void *tmp_c){

    /* Calculate the distance between two vectors given a distance matrix. */

    int i;
    int d;
    config *c;
    c = (config *)(tmp_c);

    d = 0;
    for (i = 0; i < c->num_dimensions; i++){
        d += c->dist_matrix[v1[i]][v2[i]];
    }

    return d;
}

int dl_distance(long *v1, long *v2, void *c_tmp){

    int i, j, k, max_dist;
    int db, len, minimum, insertion, deletion, transpose, cost; 
    config *c; 

    c = (config *)(c_tmp);

    max_dist = c->num_dimensions + c->num_dimensions;

    /* Populate d */
    c->dl_d[0][0] = max_dist;
    for (i = 0; i < c->num_dimensions + 1; i++){
        c->dl_d[i+1][  0] = max_dist;
        c->dl_d[i+1][  1] = i;  
    
        c->dl_d[  0][i+1] = max_dist;
        c->dl_d[  1][i+1] = i;  
    }   

    /* Calculate required moves */
    for (i = 1; i < c->num_dimensions + 1; i++){
        db = 0;
        for (j = 1; j < c->num_dimensions + 1; j++){

            k = c->dl_da[v2[j-1]];
            len = db;
            if (v1[i-1] == v2[j-1]){
                cost = 0;
                db = j;
            } else {
                cost = 1;
            }

            minimum   = c->dl_d[i  ][j  ] + cost;
            insertion = c->dl_d[i+1][j  ] + 1;
            deletion  = c->dl_d[i  ][j+1] + 1;
            transpose = c->dl_d[k  ][len] + (i - k - 1) + 1 + (j - len - 1);

            if (insertion < minimum){
                minimum = insertion;
            }

            if (deletion < minimum){
                minimum = deletion;
            }

            if (transpose < minimum){
                minimum = transpose;
            }

            c->dl_d[i+1][j+1] = minimum;

        }
        c->dl_da[v1[i-1]] = i;
    }

    /* Return final distance */
    return c->dl_d[c->num_dimensions+1][c->num_dimensions+1];

}

int query_region(int point, int *tmp_neighbors, config *c) {

    /* Look for the neighbors of a point.  Put them into tmp_neighbors and return
     * the number of neighbors found. */

    int i, neighbor_counter;
    
    neighbor_counter = 0;
    for (i = 0; i < c->num_points; i++) {

        if (i != point) {
            if (c->dist_function(c->all_points[point],c->all_points[i],(void *)c) <= c->epsilon_cutoff) {
                tmp_neighbors[neighbor_counter] = i;
                neighbor_counter++;
            }
        }
    }
    
    return neighbor_counter;
}


int local_dbscan(config *c) {

    /* Run the dbscan clustering algorithm */

    int current_cluster, i, j, k;
    int num_neighbors, num_new_neighbors;

    /* ----- Initialize arrays that will hold neighbors of points. ----- */

    int *point_neighbors, *neighbor_neighbors;

    point_neighbors = malloc(c->num_points*sizeof(int));
    if (! point_neighbors){
        fprintf(stderr,"Memory allocation error\n");
        return 1;
    }

    neighbor_neighbors = malloc(c->num_points*sizeof(int));
    if (! neighbor_neighbors){
        fprintf(stderr,"Memory allocation error\n");
        return 1;
    }
    
    for (i = 0; i < c->num_points; i++){
        point_neighbors[i] = -1;
        neighbor_neighbors[i] = -1;
    }

    /* ------ Run dbscan algorithm -------------------------------------- */

    /* Go through every point */
    current_cluster = 1;
    for (i = 0; i < c->num_points; i++) {

        /* If a point has not already been assigned to a cluster */
        if ( c->cluster_assignments[i] == -1) {

            num_neighbors = query_region(i,point_neighbors,c);
            if (num_neighbors < c->min_neighbors){
            
                /* If its neighborhood is too sparse, assign it to noise cluster */
                c->cluster_assignments[i] = 0;

            } else { 
            
                /* If not, assign it to the current cluster, then expand the
                 * cluster */
                c->cluster_assignments[i] = current_cluster;
   
                /* ---------- start expand_cluster block --------------- */ 
                num_new_neighbors = 0;
                for (j = 0; j < num_neighbors; j++) {

                    /* Go through the recently identified neighbors.  If the 
                     * neighbor isn't already in a cluster... */
                    if (c->cluster_assignments[point_neighbors[j]] == -1){

                        /* Assign it to the current cluster */
                        c->cluster_assignments[point_neighbors[j]] = current_cluster;

                        /* Look for its neighbors */
                        num_new_neighbors = query_region(point_neighbors[j],
                                                         neighbor_neighbors,c);

                        /* If these neighbors have enough neighbors, they're in 
                         * the cluster too. */
                        if (num_new_neighbors >= c->min_neighbors) {
                            for (k = 0; k < num_new_neighbors; k++){
                                if (c->cluster_assignments[neighbor_neighbors[k]] == -1){
                                    c->cluster_assignments[neighbor_neighbors[k]] = current_cluster;
                                }
                            }
                        }
                    }
                   
                }

                /* ---------- end expand_cluster block --------------- */ 

                current_cluster++;

            }
        }
    }

    /* Clean up */
    free(point_neighbors); point_neighbors = NULL;
    free(neighbor_neighbors); neighbor_neighbors = NULL;

    return 0;
}

int dbscan(long **all_points,
           long **dist_matrix,
           long *cluster_assignments, 
           int num_points,
           int num_dimensions,
           int alphabet_size,
           int epsilon_cutoff,
           int min_neighbors,
           int dist_function){

    int i, status;
    config c;

    /* Construct config data type that is used internally to keep track of points
     * and cluster assignments */ 
    c.all_points = all_points;
    c.dist_matrix = dist_matrix;
    c.cluster_assignments = cluster_assignments;

    c.num_points = num_points;
    c.num_dimensions = num_dimensions;
    c.alphabet_size = alphabet_size;
    c.epsilon_cutoff = epsilon_cutoff;
    c.min_neighbors = min_neighbors;

    /* Initially assign points to -1 (unassigned) cluster */
    for (i = 0; i < c.num_points; i++){
        c.cluster_assignments[i] = -1;
    } 
   
    if (dist_function == 0){
        c.dist_function = &simple_distance;
    } else if (dist_function == 1){
        c.dist_function = &dl_distance;
    
        /* Create array da */
        c.dl_da = (int *)malloc(c.alphabet_size*sizeof(int));
        for (i = 0; i < c.alphabet_size; i++){
            c.dl_da[i] = 0;
        }
    
        /* Create array d */
        c.dl_d = (int **)malloc((c.num_dimensions + 2)*sizeof(int *));
        for (i = 0; i < c.num_dimensions + 2; i++){
            c.dl_d[i] = (int *)malloc((c.num_dimensions + 2)*sizeof(int));
        }   

    } else {
        fprintf(stderr,"distance function not recognized\n");
        return 1;
    }
 
    /* Run dbscan */ 
    status = local_dbscan(&c);
    if (status != 0){
        fprintf(stderr,"dbscan failed.\n");
        return 1;
    }

    /* Clean up damerau levenshtein */
    if (dist_function == 1){

        for (i = 0; i < c.num_dimensions; i++){
            free(c.dl_d[i]);
        }
        free(c.dl_d);
        free(c.dl_da);
    }

    return 0;

}

/*
int main(int argc, char *argv[])
{
 
    int i, status;
   
    config c;

    c.num_points = atoi(argv[2]); 
    c.num_dimensions = atoi(argv[3]); 
    c.min_neighbors = atoi(argv[4]); 
    c.epsilon_cutoff = atoi(argv[5]);
    c.alphabet_size = 20;

    loadData(argv[1], &c);

    status = dbscan(c.all_points,
                    c.dist_matrix,
                    c.cluster_assignments, 
                    c.num_points,
                    c.num_dimensions,
                    c.alphabet_size,
                    c.epsilon_cutoff,
                    c.min_neighbors);
    
} */
