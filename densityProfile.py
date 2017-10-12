from __future__ import print_function
import mdtraj as md

# Chargement
traj = md.load('C:\Users\cyprien\Desktop\M2\Modelisation\p3_p4_p5_p8\start.pdb')
topology = traj.topology


def numberOfMoleculesPerSlice(slices):
    residues = {}
    for residue in topology.residues:
        id = residue.name + str(residue.resSeq)
        if id in residues:
            id += str("bis")
        residues[id] = []
        for atom in residue.atoms:
            # xyz[frame, atome, dimension]
            z = traj.xyz[0, atom.index, 2]
            residues[id].append(z)
    return(residues)

# for key, value in residues.items():
#     print(key)
#     print(value)

residues = numberOfMoleculesPerSlice()
print(residues['DOP1'])
print(len(residues))
print(topology)