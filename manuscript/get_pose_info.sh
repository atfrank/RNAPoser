#!/bin/bash

# script to grab poses information


rnas="1EVV 1F27 1FUF 1J7T 1LC4 1MWL 1NTA 1NTB 1O9M 1TN1 1TN2 1U8D 1YLS 1YRJ 1ZZ5 2BE0 2BEE 2CKY 2EEU 2EEV 2EEW 2ET3 2ET4 2ET5 2ET8 2F4S 2F4T 2F4U 2FCX 2FCY 2FCZ 2FD0 2G5Q 2GDI 2HOJ 2HOM 2HOO 2O3V 2O3W 2O3X 2O3Y 2OE5 2OE8 2QWY 3B4B 3B4C 3C44 3D2G 3D2V 3D2X 3DVV 3E5E 3E5F 3F2Q 3F2T 3F4H 3GX2 3GX3 3GX5 3GX6 3GX7 3IQN 3IQR 3NPQ 3SLQ 3SUH 3TD1 3TZR 3WRU 4B5R 4F8U 4F8V 4FAW 4GPW 4GPX 4GPY 4K32 4L81 4LVX 4LVY 4LW0 4P20 4P3S 4PDQ 4QLM 4QLN 4TS2 4TZX 4TZY 4WCR 4YAZ 4ZNP 5BTP 5BWS 5BXK 5C45 5KX9 6BFB"
for rna in $rnas
do
    # get info
    file="../decoys_set/poses/${rna}/rmsd.txt"
    ndecoys=`awk '{print $2}' ${file} | sort -n | tail -1`
    min=`awk '{print $3}' ${file} | sort -n | head -1 | awk '{printf "%4.2f\n", $1}'`
    max=`awk '{print $3}' ${file} | sort -n | tail -1 | awk '{printf "%4.2f\n", $1}'`
    mean=`awk '{x+=$3; next} END{print x/NR}' ${file} | awk '{printf "%4.2f\n", $1}'`
    
    echo "LOO RNA: ${rna}"
    echo "Poses: $ndecoys"
    echo "Min. RMSD: ${min} Å"
    echo "Max. RMSD: ${max} Å"
    echo "Mean RMSD: ${mean} Å"
    echo ""
done

rnas="2G5K 2O3W 2XNW 3FU2 3MUM 3Q50 3SD3 3SLM 4ERJ 4FE5 4JF2 4LX5 4NYA 4XWF 4YB0 5C7W 2B57 1F1T 2YDH 3NPN 4AOB 4KQY 4L81 5KPY"
for rna in $rnas
do
    # get info
    file="../decoys_set/poses/${rna}/rmsd.txt"
    ndecoys=`awk '{print $2}' ${file} | sort -n | tail -1`
    min=`awk '{print $3}' ${file} | sort -n | head -1 | awk '{printf "%4.2f\n", $1}'`
    max=`awk '{print $3}' ${file} | sort -n | tail -1 | awk '{printf "%4.2f\n", $1}'`
    mean=`awk '{x+=$3; next} END{print x/NR}' ${file} | awk '{printf "%4.2f\n", $1}'`
    
    echo "Validation RNA: ${rna}"
    echo "Poses: $ndecoys"
    echo "Min. RMSD: ${min} Å"
    echo "Max. RMSD: ${max} Å"
    echo "Mean RMSD: ${mean} Å"
    echo ""
done



rnas="4XWF 3SD3 4ERJ 4LX5 2XNW 1F1T 3Q50 2B57 4FE5 3MUM"
for rna in $rnas
do
    # get info
    file="../decoys_set/poses/${rna}/rmsd.txt"
    ndecoys=`awk '{print $2}' ${file} | sort -n | tail -1`
    min=`awk '{print $3}' ${file} | sort -n | head -1 | awk '{printf "%4.2f\n", $1}'`
    max=`awk '{print $3}' ${file} | sort -n | tail -1 | awk '{printf "%4.2f\n", $1}'`
    mean=`awk '{x+=$3; next} END{print x/NR}' ${file} | awk '{printf "%4.2f\n", $1}'`
    
    echo "Validation Best RNA: ${rna}"
    echo "Poses: $ndecoys"
    echo "Min. RMSD: ${min} Å"
    echo "Max. RMSD: ${max} Å"
    echo "Mean RMSD: ${mean} Å"
    echo ""
done