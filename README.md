# fast_dbscan
A lightweight, fast dbscan implementation for use on peptide strings.  It uses
pure C for the distance calculations and clustering.  This code is then wrapped
in python.

*Note*: as implemented, the software assumes all sequences have the same length.

### Installation

#### pip
```
pip3 install fast_dbscan
```

#### Development version
```
git clone https://github.com/harmslab/fast_dbscan
cd fast_dbscan
sudo python3 setup.py install
```

### Usage

#### Stand-alone
This will install a convenience program called `fast_dbscan` in the path.  This
can be invoked on the command line:

```
fast_dbscan filename
```

where `filename` is a file that contains sequences of identical length.  For 
full options, invoke:

```
fast_dbscan --help
```

#### As library
```
import fast_dbscan

d = fast_dbscan.DBScanWrapper(alphabet='amino',distance_function='dl')
d.read_file(file_with_sequences)
d.run(epsilon=1,min_neighbors=12)

# Dictionary keying cluster id to sequences
clusters = d.results
```

### Distance functions

+ `simple`: add up entries in a distance matrix based on the identies of letters
  at each column in the alignment.  Currently, the software uses hamming
  distance.  This could be easily modified to use other matricies, provided 
  distances can be calculated as integers.  The matrix is populated in 
  `DBScanWrapper.__init__`.
+ `dl`: Damerau-Levenshtein distance, allowing deletion, insertion, substitution,
  and transposition. 

### Other parameters
+ `epsilon`: the maximum distance between two samples for them to be considered
  within the same neighborhood.
+ `min_neighbors`: the minimum number of sequence neighbors required to define
  a cluster

