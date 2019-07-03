import mdtraj as md

path = '../tests/aaron_tests/complex/'

fnlist = [path + 'complex_'+str(i)+'.pdb' for i in range(1,21)]
pdb = md.load(fnlist)
pdb.save_dcd(path + 'complexes_mdtraj.dcd')
