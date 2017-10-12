from __future__ import print_function
import mdtraj as md

# Chargement
traj = md.load('C:\Users\cyprien\Desktop\M2\Modelisation\p3_p4_p5_p8\start.pdb')
topology = traj.topology

listZ = traj.xyz[0, : , 2]
minZ = min(listZ)
maxZ = max(listZ)

print(minZ)
print(maxZ)

def numberOfMoleculesPerSlice():

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

def makeSlices(nbSlices):
    zRange = maxZ/nbSlices
    zList = [minZ]
    z = zRange

    while z < maxZ:
        z += zRange
        zList.append(z)
    zList.append(maxZ)

    return(zList)


# for key, value in residues.items():
#     print(key)
#     print(value)

# residues = numberOfMoleculesPerSlice()
# print(residues['DOP1'])
# print(len(residues))

z = makeSlices(10)
print(len(z))