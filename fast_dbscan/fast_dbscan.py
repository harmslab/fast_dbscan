#!/usr/bin/env python3
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__usage__ = "fast_dbscan file_with_seqs epsilon"
__date__ = "2017-06-29"

import dbscan
import numpy as np
import sys, argparse, string

class DbscanWrapper:
    """
    Run dbscan without allocating memory for a distance matrix.  Distance is
    currently hamming distance.
    """

    def __init__(self,alphabet="amino",dist_function="hamming"):
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
        elif self.alphabet == "nucleotide":
            self._alphabet_string = "*ATCGU"
        elif self.alphabet == "latin":
            self._alphabet_string = string.ascii_letters
        else:
            raise ValueError("alphabet not recongized.")
           
        if self.dist_function == "hamming":
            self._dist_function_internal = 0
        elif self.dist_function == "dl":
            self._dist_function_internal = 1
        elif self.dist_function == "custom":
            self._dist_function_internal = 2
        else:
            err = "dist_function not recognized.\n"
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

        sequences = [l.strip() for l in lines if l.strip() != ""]

        self.load_sequences(sequences)
        

    def load_sequences(self,list_of_sequences):
        """
        Load in a collection of sequences and convert to internal integer representation.
        """ 

        if len(set([len(s) for s in list_of_sequences])) != 1:
            err = "All sequences must be the same length.\n"
            raise ValueError(err)

        self.num_points = len(list_of_sequences)
        self.num_dimensions = len(list_of_sequences[0])
        self.all_points = np.ascontiguousarray(np.zeros((self.num_points,self.num_dimensions),dtype=int))       
        for i, seq in enumerate(list_of_sequences):
            self.all_points[i,:] = np.array([self._alphabet_dict[s] for s in seq])

        self.sequences = list_of_sequences[:]
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
        argv = sys.argv[:]

    parser = argparse.ArgumentParser(description="perform dbscan clustering on a set of sequences")

    parser.add_argument("seq_file",nargs=2,help="file with sequence strings (one per line)")

    parser.add_argument("--epsilon",
                        help="neighborhood cutoff. Only sequences that differ by <= epsilon are put in a cluster",
                        type=int,
                        default=1)

    parser.add_argument("--minsize",
                        help="minimum cluster size. Clusters smaller than this are not kept. If -1, use seq_length + 1",
                        type=int,
                        default=-1)
                       
    parser.add_argument("--dist",
                        help="distance function. hamming (default) or dl (Damerau-Levenshtein).  dl distance is *much* slower than hamming.",
                        type=str,
                        choices=['hamming','dl'],
                        default="hamming")

    parser.add_argument("--alphabet",
                        help="alphabet for sequence.  amino (protein seq), nucleotide (nucleotide seq), latin (upper and lowercase latin)",
                        type=str,
                        choices=['amino','nucleotide','latin'],
                        default="amino")
 
    args = parser.parse_args(argv)

    seq_file = args.seq_file[1]
    print(seq_file)
    epsilon = args.epsilon
    if epsilon < 1:
        err = "epsilon must be greater than zero.\n"
        raise ValueError(err)

    min_size = args.minsize
    if min_size == -1:
        min_size = None
    elif min_size < 2:
        err = "minsize must be > 2 OR -1.  If -1, minimum cluster size will be seq_length + 1\n"
        raise ValueError(err)
    else:
        pass
    
    dist = args.dist 
   
    alphabet = args.alphabet 

    d = DbscanWrapper(alphabet=alphabet,dist_function=dist)
    d.read_file(seq_file)
    d.run(epsilon=epsilon,min_neighbors=min_size)

    clusters = d.results
    cluster_ids = list(clusters.keys())
    for c in cluster_ids:
        seqs = clusters[c]
        for s in seqs:
            print(c,s)

if __name__ == "__main__":
    main()
