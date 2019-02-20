/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Jingru Xie and Aaron T. Frank
     
*/

#include "Molecule.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"
#include "LARMORD.hpp"
#include "Trajectory.hpp"
#include "AtomicFeaturizer.hpp"
#include "MolecularFeaturizer.hpp"
#include "../lib/Metallo.hpp"


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <time.h> // keep track of processing time of program

using namespace std;

void usage(){
  std::cerr << "====================================================" << std::endl;
  std::cerr << "ATOMIC Featurizer v1.00" << std::endl;
  std::cerr << "(c) 2017 Jingru Xie, Aaron T. Frank and University of Michigan." << std::endl;
  std::cerr << "====================================================" << std::endl;
  std::cerr << "Usage:   rna_poser [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-mode prediction mode 'R','L' or 'RL']" << std::endl;
  std::cerr << "         [-outfile path and name (without extension) of output file]" << std::endl;
  std::cerr << "         [-mol2 MOL2file]" << std::endl;
  std::cerr << "         [-rdock rdock score file]" << std::endl;
  std::cerr << "         [-trj TRAJfile]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;  
  std::cerr << "         [-identification ID]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){
  int i;
  unsigned int f;
  clock_t t_start, t_end;

  t_start = clock();
  
  std::stringstream resid;
  std::vector<std::string> pdbs;
  std::vector<std::string> mol2s;
  std::vector<std::string> rdocks;
  std::string currArg;
  std::string nucleus;
  std::string resname;
  std::string resnameCode;
  std::string atomname;
  std::string key, residID;
  std::string outf;
  std::string outfile;
  std::vector<int> selected_residues;
  std::vector<std::string> selected_nuclei;
  bool isResidue;
  bool mismatchCheck;
  std::map<std::string, double> histo;
  std::string moltype;
  
  std::vector<std::string> trajs;
  int start;
  int stop=std::numeric_limits<int>::max();
  int skip;
  bool startFlag=false;
  unsigned long long int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned long long int nframe;
  unsigned long long int process;

  start=0;
  skip=0;
  nframe=0;
  process=0;
  
  int beta;
  
  bool scalar = false;
  bool molecular = true;
  bool normalization = true;
  double cutoff;
  double eta;
  int numEta;
  int etaStartPow;
  double etaBase;
  string select_atm;
  vector<string> sel_atmname;
  vector<string> row_atmname;
  size_t pos = 0;
  string atm;
  string delim = " ";
  // hard code sel_atmname
select_atm=":ADE.C1' :ADE.C2 :ADE.C2' :ADE.C3' :ADE.C4 :ADE.C4' :ADE.C5 :ADE.C5' :ADE.C6 :ADE.C8 :ADE.N1 :ADE.N3 :ADE.N6 :ADE.N7 :ADE.N9 :ADE.O2' :ADE.O3' :ADE.O4' :ADE.O5' :ADE.OP1 :ADE.OP2 :ADE.P :CYT.C1' :CYT.C2 :CYT.C2' :CYT.C3' :CYT.C4 :CYT.C4' :CYT.C5 :CYT.C5' :CYT.C6 :CYT.N1 :CYT.N3 :CYT.N4 :CYT.O2 :CYT.O2' :CYT.O3' :CYT.O4' :CYT.O5' :CYT.OP1 :CYT.OP2 :CYT.P :GUA.C1' :GUA.C2 :GUA.C2' :GUA.C3' :GUA.C4 :GUA.C4' :GUA.C5 :GUA.C5' :GUA.C6 :GUA.C8 :GUA.N1 :GUA.N2 :GUA.N3 :GUA.N7 :GUA.N9 :GUA.O2' :GUA.O3' :GUA.O4' :GUA.O5' :GUA.O6 :GUA.OP1 :GUA.OP2 :GUA.P :URA.C1' :URA.C2 :URA.C2' :URA.C3' :URA.C4 :URA.C4' :URA.C5 :URA.C5' :URA.C6 :URA.N1 :URA.N3 :URA.O2 :URA.O2' :URA.O3' :URA.O4 :URA.O4' :URA.O5' :URA.OP1 :URA.OP2 :URA.P";
  while ((pos = select_atm.find(delim)) != string::npos){
    atm = select_atm.substr(0, pos);
    sel_atmname.push_back(atm);
    select_atm.erase(0, pos + delim.length());
  }
  sel_atmname.push_back(select_atm);
  // hard code row_atmname
  row_atmname.push_back(":UNK.");

  cutoff = 20.0;
  numEta = 3;
  etaStartPow = 1;
  etaBase = 2.0;
  string mode="R";
 
  Molecule *neighbormol;
  neighbormol=NULL;
  
  Atom *ai, *aj;  
  ai=NULL;
  aj=NULL;
  outfile="prediction.txt";
  beta=-3;
  moltype = "protein";
  
  Metallo *metallo = new Metallo(1);
  vector<double> feature_tmp;
  vector<int> pred;
  vector<double> negprob;
  vector<double> prob;
  int numrows = 0;
  double rmsdcutoff = 2.5;
  double probability;
  double negprobability;
  int prediction;
 
  LARMORD *larm;
  larm=NULL;

  pdbs.clear();
  selected_residues.clear();
  selected_nuclei.clear();
  isResidue = true;
  mismatchCheck = false;
  histo.clear();
  
  for (int i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0)
    {
      usage();
    }
	else if (currArg.compare("-mode") == 0 ){
      currArg = argv[++i];
      mode = currArg;
    }
    else if (currArg.compare("-rmsdcutoff") == 0 ){
      currArg = argv[++i];
      std::stringstream(currArg)>>rmsdcutoff;
    }
    else if (currArg.compare("-rowatm") == 0){
      currArg=argv[++i];
      select_atm = currArg;
      size_t pos = 0;
      string atm;
      string delim = " ";
      while ((pos = select_atm.find(delim)) != string::npos){
        atm = select_atm.substr(0, pos);
        row_atmname.push_back(atm);
        select_atm.erase(0, pos + delim.length());
      }
      row_atmname.push_back(select_atm);
    }
    else if (currArg.compare("-outfile") == 0)
    {
      currArg=argv[++i];
      outfile = currArg;
    }
    else if (currArg.compare("-mol2") == 0)
    {
      currArg=argv[++i];
      mol2s.push_back(currArg);
    }
    else if (currArg.compare("-rdock") == 0)
    {
      currArg=argv[++i];
      rdocks.push_back(currArg);
    }
    else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0)
    {
      currArg=argv[++i];
      trajs.push_back(currArg);
    }
    else if (currArg.compare("-skip") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
      startFlag=true;
    }
    else if (currArg.compare("-stop") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }    
    else if (currArg.compare(0,1,"-") == 0)
    {
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }
  if (pdbs.size() == 0 && mode != "L")
  {
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }
  if (rdocks.size() == 0 && mode != "R")
  {
  std::cerr << std::endl << "Error: Please provide a rdock score file" << std::endl << std::endl;
  usage();
  }
  
  // Case L: predict based on rdock score only
  if (mode == "L"){
    std::string FILENAME = rdocks.at(0);
    std::ifstream infile(FILENAME);
    std::string line;
    double temp;
    int prediction;
    numrows = 0;

    while (std::getline(infile, line))
    {
      numrows = numrows + 1;
      std::istringstream iss(line);
      feature_tmp.clear();
      // iss >> temp;
      for(i = 0; i < 30; i++){
        iss >> temp;
        feature_tmp.push_back(temp);
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
    }
    ofstream outFile;
    outFile.open(outfile.c_str());
    cout << "==============" << endl;
    cout << "pred p0 p1" << endl;
    outFile << "pred p0 p1" << endl;
    for(i = 0; i < numrows; i ++){
      cout << pred.at(i) << " " << negprob.at(i) << " " << prob.at(i) << endl;
      outFile << pred.at(i) << " " << negprob.at(i) << " " << prob.at(i) << endl;
    }
    cout << endl;
    return 0;
  }
  
  //initialize
  Molecule *mol=NULL;
  Molecule *mol2=NULL;
  std::vector<Atom*> atmVec;
  AtomicFeaturizer *atomicfeature;
  atomicfeature=NULL;
  std::vector<string> atmType;
  
	std::vector<double> etalist;
	std::vector<std::vector<double> > features;
  std::vector<double> mfeats;
  std::vector<std::vector<double> > mol_features;
  std::vector<int> frames; // frameids being featurized in trajectory analysis
  
  // initialize SYBYL sequence
  std::string str[] = {"C.1","C.2","C.3","C.ar","C.cat","H","N.1","N.2","N.3","N.4","N.am","N.ar","N.pl3","O.2","O.3","O.co2","P.3","S.2","S.3","S.o","S.o2"};
  const std::vector<string> SYBYL(str, str+(sizeof(str)/sizeof(string)));
  
  //create list of eta values
	for (int i = 0; i < numEta; i ++){
		eta = pow(etaBase, i + etaStartPow);
		etalist.push_back(eta);
	}

  
  if (trajs.size() > 0){
    if (pdbs.size() > 1){
      std::cerr << std::endl << "Warning: Only the first PDB structure is used for trajectory analysis" << std::endl << std::endl;
      }
    /* Trajectory analysis */
    mol=Molecule::readPDB(pdbs.at(0));
    mol->selAll();
    
    // initialize mol2 atmtypes for molecular featurization
    if (molecular){
      if (!mol2s.size()){
        cout << endl <<  "Error: Please provide a mol2 file. [-mol2 MOL2file] " << endl << endl;
        return 1;
      }
      mol2 = getAtmType(mol2s.at(0), atmType);
    }
    
    /* Process trajectories */
    for (itrj=0; itrj< trajs.size(); itrj++){
      trjin.open(trajs.at(itrj).c_str(), std::ios::binary);
      frames.clear();
      mol_features.clear();
      if (trjin.is_open()){
        ftrjin=new Trajectory;
        ftrjin->setMolecule(mol);
        if (ftrjin->findFormat(trjin) == true){
          ftrjin->readHeader(trjin);
          if (skip > 0 && startFlag == false){
            start=skip;
            }
          /* Print out current molecule info */
          cout << trajs.at(itrj) << endl;
          cout << "Number of Atoms: " << mol->getAtmVecSize() << endl;
          cout << "Number of Residues: " << mol->getResVecSize() << endl;
          cout << "Number of Chains: " << mol->getChnVecSize() << endl;
          cout << "Number of frames in traj: " << ftrjin->getNFrame() << endl;
          /* Loop through desired frames */
          for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip){
            if (ftrjin->readFrame(trjin, i) == false){ // read frame and setmol here
              std::cerr << "Warning: EOF found before the next frame could be read" << std::endl;
              break;
              }
            nframe++;
            frames.push_back(i + 1);
            cout << endl << "Frame " << i+1 << ":" << endl;
            atomicfeature = new AtomicFeaturizer(mol);
            if (molecular){
              atomicfeature->featurizeScalar(cutoff, etalist, features, outf, sel_atmname, row_atmname, false);
              MolecularFeaturizer(mol, mol2, features, mol_features, outf, atmType, SYBYL, mfeats, false);
            }
            else if (scalar){
              atomicfeature->featurizeScalar(cutoff, etalist, features, outf, sel_atmname, row_atmname);
            }
            else{
              atomicfeature->featurize(cutoff, etalist, features,
                                       outf, sel_atmname, row_atmname);
            }
            features.clear();
          }
        }
        else{
          std::cerr << "Warning: Skipping unknown trajectory format \"";
          std::cerr << trajs.at(itrj) << "\"" << std::endl;
        }
        if (ftrjin != NULL){
          delete ftrjin;
        }
        if (atomicfeature != NULL){
          delete atomicfeature;
        }
        if (molecular){
          if (normalization){
            normalize(mol_features);
          }
        }
      }
      trjin.close();
    }
  }

  else{
    if (molecular){
       /* molecular featurization */
      if (!mol2s.size()){
        cout << endl <<  "Error: Please provide a mol2 file. [-mol2 MOL2file] " << endl << endl;
        return 1;
      }
      mol2 = getAtmType(mol2s.at(0), atmType);
    }
    for (f=0; f< pdbs.size(); f++){
      mol=Molecule::readPDB(pdbs.at(f));
      /* Print out file and atom info */
      atomicfeature = new AtomicFeaturizer(mol);
      cout << "Input: " << pdbs.at(f) << endl;
      cout << "Number of Atoms: " << mol->getAtmVecSize() << endl;
      cout << "Number of Residues: " << mol->getResVecSize() << endl;
      cout << "Number of Chains: " << mol->getChnVecSize() << endl;
      outf = outfile + ".txt";
      
      if (molecular){
        if (mol2s.size() > 1){
          mol2 = getAtmType(mol2s.at(f), atmType);
        }
        atomicfeature->featurizeScalar(cutoff, etalist, features, outf, sel_atmname, row_atmname, false);
        MolecularFeaturizer(mol, mol2, features, mol_features, outf, atmType, SYBYL, mfeats, false);
      }
      else if (scalar){
        atomicfeature->featurizeScalar(cutoff, etalist, features, outf, sel_atmname, row_atmname);
      }
      else{
        atomicfeature->featurize(cutoff, etalist, features,
                                   outf, sel_atmname, row_atmname);
      }
      features.clear();
    }
    if(atomicfeature!=NULL){
      delete atomicfeature;
    }
  }
  
  t_end = clock();
  float diff ((float)t_end - (float)t_start);
  cout << "Processing time: " << diff/CLOCKS_PER_SEC << " seconds." << endl << endl;

  numrows = 0;
  // Case RL: predict based on mol features + rdock score
  if (mode == "RL"){
    std::string FILENAME = rdocks.at(0);
    std::ifstream infile(FILENAME);
    std::string line;
    double temp;
    int prediction;
    numrows = 0;
    int j = 0;
    
    while (std::getline(infile, line))
    {
      numrows = numrows + 1;
      if (numrows < start || numrows < skip || (numrows-start) % (skip+1) != 0){
        continue;
      }
      if (numrows > stop){
        break;
      }
      std::istringstream iss(line);
      feature_tmp.clear();
      // mol features
      feature_tmp = mol_features.at(j);
      j = j + 1;
      // rdock scores
      // iss >> temp;
      for(i = 0; i < 30; i++){
        iss >> temp;
        feature_tmp.push_back(temp);
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
    }
    ofstream outFile;
    outFile.open(outfile.c_str());
    cout << "==============" << endl;
    cout << "pred p0 p1" << endl;
    outFile << "pred p0 p1" << endl;
    for(i = 0; i < j; i ++){
      cout << pred.at(i) << " " << negprob.at(i) << " " << prob.at(i) << endl;
      outFile << pred.at(i) << " " << negprob.at(i) << " " << prob.at(i) << endl;
    }
    cout << endl;
    return 0;
  }

  numrows = 0;
  // Case R: predict based on mol features only
  for (i = 0; i < mol_features.size(); i ++)
  {
	  feature_tmp = mol_features.at(i);
    numrows = numrows + 1;
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
  }
  ofstream outFile;
  outFile.open(outfile.c_str());
  cout << "==============" << endl;
  cout << "pred p0 p1" << endl;
  outFile << "pred p0 p1" << endl; 
  for(i = 0; i < numrows; i ++){
    cout << pred.at(i) << " " << negprob.at(i) << " " << prob.at(i) << endl;    
    outFile << pred.at(i) << " " << negprob.at(i) << " " << prob.at(i) << endl;    
  }
  cout << endl;
  return 0;
};
