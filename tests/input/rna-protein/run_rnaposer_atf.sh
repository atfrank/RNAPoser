#!/bin/bash

# RNAPOSER
module load gcc
module load anaconda
base=/home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/
rnaposer=/home/itssahil/PROJECTS/Mol_Feturizer/Test_systems_manuscript/RNAPoser-master/bin/rna_poser
pdbs="1B7F"
for pdb in $pdbs
do
    ${rnaposer} ${base}/${pdb}_complex.pdb -trj ${base}/${pdb}_complex.dcd -mol2 ${base}/${pdb}_protein_from_pdb.mol2 -mode R -outfile ${base}/${pdb}_classifications.txt
done    


# Python
./src/rna_poser_validation.sh ${base}/ ${pdb} ${pdb}_complex.pdb ${pdb}_protein_from_pdb.mol2 ${pdb}_complex.dcd
