#!/bin/bash
if [ "$#" -ge 2 ] && [ "$1" = "R" ]
then
	cp predictors_RF/predictors_${2}/R_all.cpp lib/Metallo_RF.cpp
else
	cp predictors_RF/${1}_all.cpp lib/Metallo_RF.cpp
fi
make

