#!/bin/bash

ntest=1
cd tests
test=1

../bin/featurize -cutoff 5.0 -eta 1.0 -numEta 3 -etaStartPow 0 -etaBase 2.0 -output cluster02_1.dat.test cluster02_1.pdb  > extractor_test.log
paste <(cat python_cluster02_1_comparison.dat) <(cat cluster02_1.dat.test) | awk -v test=$test '{print ($1-$4)+($2-$5)+($3-$6)}' | awk '{ sum += $1; n++} END {print sum}' | awk '{ if($1<0.0001) print "TEST PASS "test; else print "TEST FAIL "test}'

