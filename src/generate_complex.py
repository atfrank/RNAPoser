import sys
import inspect
from glob import glob
from pymol import cmd
import os

def split_decoys(poses):
    for a in range(1, 1+cmd.count_states("(%s)"%poses)):
        cmd.frame(a)
        cmd.save("poses/poses_%i.pdb"%a)
        if a == 1:
            cmd.save("lig.sd")

def generate_complexes(receptor="receptor.mol2", poses="poses.sd"):
    os.system("mkdir -p complex")
    os.system("mkdir -p poses")
    cmd.delete("*")
    cmd.load(poses, "poses")
    split_decoys("poses")
    poses=range(1, 1+cmd.count_states("(poses)"))
    cmd.delete("poses")
    for a in poses:
        cmd.delete("*")
        cmd.load(receptor)
        cmd.load("poses/poses_%i.pdb"%a)
        cmd.create("test", "all")
        cmd.save("complex/complex_%i.pdb"%a,"test")

generate_complexes(receptor="receptor.pdb", poses="poses.sd")
