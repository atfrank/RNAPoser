#!/bin/bash
#SBATCH -p frank
#SBATCH -A frank
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
if [[ $# -ne 1 ]]
then
    echo "usage: $0 <receptor index (e.g., 0, 1, 2,...)>"
else
    receptor=$1
    ligands=`seq 1 9810`
    rm -f reports/rna-poser-scores.txt
    for ligand in ${ligands}
    do
        pdbid=YALE-${receptor}-${ligand}
        if [[ -f rna_poser/${pdbid}/poser-scores.txt ]]
        then
            if [[ -f scores/docking_titled_${receptor}_${ligand}.sd ]]
            then
                 ligand_file="scores/docking_titled_${receptor}_${ligand}.sd"
            else
                 ligand_file="scores/docking_${receptor}_${ligand}.sd"
            fi
            
            ligand_name1=`head -1 ${ligand_file}`
            ligand_name2=`head -1 ${ligand_file}`
            awk -v ligand_number=${ligand} -v lig1=${ligand_name1} -v lig2=${ligand_name2} '{print ligand_number, lig1, lig2, $0}' rna_poser/${pdbid}/poser-scores.txt | tee -a reports/rna-poser-scores.txt
        fi
    done
fi