#!/bin/bash

for i in `cat rnas_C.txt`;
do
  cp ../predictors/RF_${i}.predictor.pkl predictor.pkl
  sh getTrees.bash 1000 > ../predictors/${i}_RF.cpp
done
