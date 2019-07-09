#!/bin/bash

# usage: ./src/rna_poser.sh tests/input/${pdb}/ ${pdb} receptor.pdb poses.sd
source_dir=$1

pdbid=$2
receptor="$1$3"
poses="$1$4"
rmsd=$5
eta=$6

if [[  -f "${poses}" ]]
then
    # check if rmsd valid
    if ! [[ ":1:1.5:2:2.5:" = *:$rmsd:* ]]
    then
        echo "ERROR: rmsd value not valid (need to be 1, 1.5, 2, 2.5). Exit."
        exit
    fi
    # check if eta valid
    if [ "$eta" == "2" ]
    then
        nEta=1
    elif [ "$eta" == "24" ]
    then
        nEta=2
    elif [ "$eta" == "248" ]
    then
        nEta=3
    else
        echo "ERROR: eta value not valid (need to be 2, 24 or 248). Exit."
        exit
    fi

    # make direct to store data
    rm -rf working_dir/${pdbid}/
    mkdir -p working_dir/${pdbid}/
    echo $poses

    # copy need data
    cp ${receptor} working_dir/${pdbid}/receptor.pdb
    cp ${poses} working_dir/${pdbid}/poses.sd

    # get complexes
    cd working_dir/${pdbid}
    $PYMOL -cqr ../../src/generate_complex.py #> /dev/null
    echo "finished generating complexes for ${pdbid}"

    # fix complex names
    grep ATOM complex/complex_1.pdb > tmp_rec.pdb
    grep HETATM complex/complex_1.pdb > tmp_lig.pdb
    cat tmp_lig.pdb | awk '{print "UNK Z   1  "}' > lignames.txt
    cat tmp_rec.pdb | awk '{print $4}' | \
    awk '{ if(substr($1, 0, 1)=="G") $1="GUA"; print $1}' | \
    awk '{ if(substr($1, 0, 1)=="A") $1="ADE"; print $1}' | \
    awk '{ if(substr($1, 0, 1)=="C") $1="CYT"; print $1}' | \
    awk '{ if(substr($1, 0, 1)=="U" && $1!="UNK") $1="URA"; print $1}' > resnames.txt
    paste -d "\0" <(cut -c1-17 tmp_rec.pdb) <(cat resnames.txt) <(cut -c21- tmp_rec.pdb) > complex_receptor.pdb
    paste -d "\0" <(cut -c1-17 tmp_lig.pdb) <(cat lignames.txt) <(cut -c29- tmp_lig.pdb) > complex_ligand.pdb
    cat <(cat complex_receptor.pdb) <(echo "TER") <(cat complex_ligand.pdb) <(echo "END") > complex.pdb
    rm -f tmp_lig.pdb tmp_rec.pdb lignames.txt resnames.txt complex_receptor.pdb complex_ligand.pdb

    # make dcd file
    decoys=`ls complex/complex_*.pdb | wc -l | awk '{print $1}'`
    seq 1 ${decoys} | awk '{print "complex/complex_"NR".pdb"}' > pdblist
    perl /Users/jingrux/Documents/GitSoftware/toolset-master/perl/pdb2traj.pl -f pdblist -out complexes.dcd

    # get ligand info
    cat complex.pdb | sed 's/O1P/OP1/g' | sed 's/O2P/OP2/g' | sed 's/C\.2/C2 /g' | sed 's/N\.1/N1 /g' | sed 's/N\.2/N2 /g' | sed 's/N.3/N3 /g' | sed 's/N\.4/N4 /g' | sed 's/O\.2/O2 /g' | sed 's/P\.3/P  /g' | sed 's/OP3/OP3 /g' > complex_fix.pdb

    mv complex_fix.pdb complex.pdb
    grep 'UNK' complex.pdb > lig.pdb
    babel -isd lig.sd -omol2 lig.mol2

    # featurize trajectory
    cd ../..
    selatms=":ADE.C1' :ADE.C2 :ADE.C2' :ADE.C3' :ADE.C4 :ADE.C4' :ADE.C5 :ADE.C5' :ADE.C6 :ADE.C8 :ADE.N1 :ADE.N3 :ADE.N6 :ADE.N7 :ADE.N9 :ADE.O2' :ADE.O3' :ADE.O4' :ADE.O5' :ADE.OP1 :ADE.OP2 :ADE.P :CYT.C1' :CYT.C2 :CYT.C2' :CYT.C3' :CYT.C4 :CYT.C4' :CYT.C5 :CYT.C5' :CYT.C6 :CYT.N1 :CYT.N3 :CYT.N4 :CYT.O2 :CYT.O2' :CYT.O3' :CYT.O4' :CYT.O5' :CYT.OP1 :CYT.OP2 :CYT.P :GUA.C1' :GUA.C2 :GUA.C2' :GUA.C3' :GUA.C4 :GUA.C4' :GUA.C5 :GUA.C5' :GUA.C6 :GUA.C8 :GUA.N1 :GUA.N2 :GUA.N3 :GUA.N7 :GUA.N9 :GUA.O2' :GUA.O3' :GUA.O4' :GUA.O5' :GUA.O6 :GUA.OP1 :GUA.OP2 :GUA.P :URA.C1' :URA.C2 :URA.C2' :URA.C3' :URA.C4 :URA.C4' :URA.C5 :URA.C5' :URA.C6 :URA.N1 :URA.N3 :URA.O2 :URA.O2' :URA.O3' :URA.O4 :URA.O4' :URA.O5' :URA.OP1 :URA.OP2 :URA.P"
    rowatms=":UNK."

    ./bin/featurize \
    -etaStartPow 1 \
    -numEta $nEta \
    -cutoff 20 \
    -outfile working_dir/${pdbid}/features_ \
    -scalar 1 \
    -molecular 1 \
    -mol2 working_dir/${pdbid}/lig.mol2 \
    -normalization 1 \
    -rowatm "`echo ${rowatms}`" \
    -selatm "`echo ${selatms}`"\
    -trj working_dir/${pdbid}/complexes.dcd working_dir/${pdbid}/complex.pdb

    # get score
    python src/rna_poser.py --classifier classifier/eta${eta}/RF_${rmsd}.predictor.pkl  --features working_dir/${pdbid}/features_traj1.txt --output working_dir/${pdbid}/poser-scores.txt
fi
