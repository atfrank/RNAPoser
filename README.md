# RNAPoser
Machine-Learning Pose Classifier for RNA-Ligand Complexes

## Install
```shell
$ cd /path/to/RNAPoser/
$ make clean
$ make
```

## Usage manual
```shell
$ sh make_predictor.sh [mode] [RMSD_threshold (optional, only for R mode, possible values: {10, 15, 20, 25, 30}(*0.1 A))]
$ ./bin/rna_poser -h
Usage:   rna_poser [-options] <PDBfile>
Options: [-mode prediction mode 'R','L' or 'RL']
         [-outfile path and name (without extension) of output file]
         [-mol2 MOL2file]
         [-rdock rdock score file]
         [-trj TRAJfile]
         [-skip frames] [-start frame] [-stop frame]
         [-identification ID]

```
## example
```shell
$ sh make_predictor.sh R
$ ./bin/rna_poser sahil-tests/complex.pdb -trj sahil-tests/complexes.dcd -mol2 sahil-tests/lig_2b57.mol2 -mode R

file: prediction.txt
pred p0 p1
1 0.01 0.99
1 0.002 0.998
1 0.002 0.998
1 0 1
1 0 1
1 0.001 0.999
1 0 1
1 0 1
1 0 1
1 0 1
1 0.001 0.999
  ...
```
```shell
$ sh make_predictor.sh R 10
$ ./bin/rna_poser sahil-tests/complex.pdb -trj sahil-tests/complexes.dcd -mol2 sahil-tests/lig_2b57.mol2 -mode R

file: prediction.txt
pred p0 p1
1 0.033 0.967
1 0.007 0.993
1 0.002 0.998
1 0.002 0.998
1 0.002 0.998
1 0.001 0.999
1 0.002 0.998
1 0.001 0.999
1 0.001 0.999
1 0 1
  ...
```

```shell
$ sh make_predictor.sh RL
$ ./bin/rna_poser sahil-tests/complex.pdb -trj sahil-tests/complexes.dcd -mol2 sahil-tests/lig_2b57.mol2 -mode RL -rdock sahil-tests/Scores.txt

file: prediction.txt
1 0.101 0.899
1 0.082 0.918
1 0.077 0.923
1 0.076 0.924
1 0.074 0.926
1 0.074 0.926
1 0.074 0.926
1 0.074 0.926
1 0.074 0.926
1 0.074 0.926
  ...
```


## License
```
  Copyright University of Michigan.
  Author: Jingru Xie and Aaron T. Frank

```
