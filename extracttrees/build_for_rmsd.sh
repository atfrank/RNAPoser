#!/bin/bash

for i in `cat rmsds.txt`;
do
  cp trees/RF_${i}.predictor.pkl  predictor.pkl
  sh getTrees.bash 1000 > ../predictors_RF/predictors_${i}/R_all.cpp
done
