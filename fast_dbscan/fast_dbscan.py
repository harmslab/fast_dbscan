#!/usr/bin/env python3
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__usage__ = "fast_dbscan file_with_seqs epsilon"
__date__ = "2017-06-29"

import dbscan
import numpy as np
import sys

class DbscanWrapper:
    """
    Run dbscan without allocating memory for a distance matrix.  Distance is
    currently hamming distance.
    """

    def __init__(self,alphabet="amino",dist_function="simple"):
        """
        Initialize the class.  This should be called by every subclass to
        initialize the internal dictionaries mapping alphabets to fast internal
        indexes. 
        """

        # initialize internal variables
        self.alphabet = alphabet
        self.dist_function = dist_function

        # decide on the alphabet
        if self.alphabet == "amino": 
            self._alphabet_string = "*ABCDEFGHIKLMNPQRSTVWXYZ"
        else:
            raise ValueError("alphabet not recongized.")
           
        if self.dist_function == "simple":
            self._dist_function_internal = 0
        elif self.dist_function == "dl":
            self._dist_function_internal = 1
        else:
            err = "dist_function not recognized.  should be 'simple' or 'dl' (Damerau-Levenshtein)\n"
            raise ValueError(err)
 
        self.alphabet_size = len(list(self._alphabet_string))

        enum_list = zip(self._alphabet_string,range(len(self._alphabet_string)))
        self._alphabet_dict = dict([(a, i) for a, i in enum_list])

        tmp_matrix = np.zeros((self.alphabet_size,self.alphabet_size),dtype=int)
        for k1 in self._alphabet_string:
            i = self._alphabet_dict[k1]   
            for k2 in self._alphabet_string:
                j = self._alphabet_dict[k2]
                if k1 == k2:
                    tmp_matrix[i,j] = 0
                else:
                    tmp_matrix[i,j] = 1

        self.dist_matrix = tmp_matrix

    def read_file(self,filename):
        """
        Read file with sequences and convert to integer representation.
        """

        f = open(filename,'r')
        lines = f.readlines()
        f.close()

        self.num_points = len(lines)
        self.num_dimensions = len(lines[0].strip())
        self.all_points = np.ascontiguousarray(np.zeros((self.num_points,self.num_dimensions),dtype=int))       
        for i, l in enumerate(lines):
            self.all_points[i,:] = np.array([self._alphabet_dict[s] for s in l.strip()])

        self.sequences = [l.strip() for l in lines]
        self.cluster_assignments = np.ascontiguousarray(np.ones(self.num_points,dtype=int)*-1)

    def run(self,epsilon,min_neighbors=None):
        """
        Run calculation.
        """

        if min_neighbors is None:
            min_neighbors = self.num_dimensions + 1

        status = dbscan.run_dbscan(self.all_points,
                                   self.dist_matrix,
                                   self.cluster_assignments,
                                   self.num_points,
                                   self.num_dimensions,
                                   self.alphabet_size,
                                   epsilon,
                                   min_neighbors,
                                   self._dist_function_internal)

    @property
    def results(self):

        clusters = {}
        for i in range(len(self.cluster_assignments)):
            try:
                clusters[self.cluster_assignments[i]].append(self.sequences[i])
            except KeyError:
                clusters[self.cluster_assignments[i]] = [self.sequences[i]]
                
        cluster_ids = np.unique(self.cluster_assignments)

        return clusters

def main(argv=None):
    
    if argv is None:
        argv = sys.argv[1:]

    try:
        file_name = argv[0]
        epsilon = int(argv[1])
    except (IndexError,ValueError):
        err = "Incorrect arguments. Usage:\n\n{}\n\n".format(__usage__)
        raise ValueError(err)

    try:
        dist_function = argv[2]
    except IndexError:
        dist_function = "simple"

    d = DbscanWrapper(dist_function=dist_function)
    d.read_file(file_name)
    d.run(epsilon)

    clusters = d.results
    cluster_ids = list(clusters.keys())
    for c in cluster_ids:
        seqs = clusters[c]
        for s in seqs:
            print(c,s)

if __name__ == "__main__":
    main()
