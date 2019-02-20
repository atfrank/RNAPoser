/*
 
 Copyright University of Michigan.
 This file is part of the Larmor software suite and is made available under license.
 University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.

 */

#include "../lib/Metallo.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <iomanip>

using namespace std;

void usage(){
    cerr << "====================================================" <<endl;
    cerr << "Predict" << std::endl;
    cerr << "====================================================" <<endl;
}

int main(int argc, char **argv){
  Metallo *metallo = new Metallo(1);
  vector<double> feature_tmp;
  vector<int> pred;
  vector<double> negprob;
  vector<double> prob;
  vector<int> rmsdtag;
  vector<double> rmsd;
  int numrows = 0;
  string rnaid = "1aju";
  double rmsdcutoff = 2.5;
  double temp;
  double probability;
  double negprobability;
  int prediction;
  int i;
  string currArg;
  
  // input parameters
  for (int i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-rnaid") == 0 ){
      currArg = argv[++i];
      rnaid = currArg;
    }
    else if (currArg.compare("-rmsdcutoff") == 0 ){
      currArg = argv[++i];
      std::stringstream(currArg)>>rmsdcutoff;
    }
  }
  
  std::string FILENAME = rnaid + "_test.txt";
  std::ifstream infile(FILENAME);
  std::string line;
  while (std::getline(infile, line))
  {
    numrows = numrows + 1;
    std::istringstream iss(line);
    feature_tmp.clear();
    iss >> temp;
    for(i = 0; i < 5387; i++){
      iss >> temp;
      if (i < 5355){
        feature_tmp.push_back(temp);
      }
      if (i == 5386){
        rmsd.push_back(temp);
        if (temp <= rmsdcutoff){
          rmsdtag.push_back(1);
        }
        else{
          rmsdtag.push_back(0);
        }
      }
    }
    probability = metallo->Predict(feature_tmp);
    negprobability = 1 - probability;
    prob.push_back(probability);
    negprob.push_back(negprobability);
    if (probability >= 0.5){
      prediction = 1;
    }
    else{
      prediction = 0;
    }
    pred.push_back(prediction);
  
    //cout << feature_tmp.at(0) << feature_tmp.at(1) << endl;
    // process pair (a,b)
  }
  for(i = 0; i < numrows; i ++){
    cout << pred.at(i) << " " << negprob.at(i) << " " << prob.at(i) << " " << rmsdtag.at(i) << " " << rmsd.at(i) << endl;
    
  }

return 0;
};
