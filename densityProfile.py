#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Chargement trajectoires/topologie
traj = md.load('C:\Users\cyprien\Desktop\M2\Modelisation\p3_p4_p5_p8\md_200ns_OK.xtc',
               top='C:\Users\cyprien\Desktop\M2\Modelisation\p3_p4_p5_p8\start.pdb')
topology = traj.topology

def sort_atoms_inside_slices(slices, frame_id):
    h2o = {i: [] for i in range(len(slices))}
    dopc = {i: [] for i in range(len(slices))}
    dope = {i: [] for i in range(len(slices))}

    for residue in topology.residues:
        list_pos_z_atoms = []
        is_dope = False
        for atom in residue.atoms:
            # Si atome d'H et appartient a un DOP --> DOPE
            if atom.element.symbol == 'H' and not atom.residue.is_water:
                is_dope = True
            # Recuperation position z de l'atome (xyz[frame, indice de l'atome, dimension])
            list_pos_z_atoms.append(traj.xyz[frame_id, atom.index, 2])
        # Parcours des tranches pour placer les atomes dans les bons intervals
        # Regle: indice 1 (tranche 1) de valeur z1 --> contient tous les atomes
        # de la tranche 1 qui auront un z compris entre z0 (le min) et z1.
        for z in list_pos_z_atoms:
            i = 1
            while z > slices[i] and i < len(slices) - 1:
                i += 1
            # Atome appartenant à une molecule H2O
            if(atom.residue.is_water):
                h2o[i].append(z)
            # Atome appartenant a un DOPE
            elif(is_dope):
                dope[i].append(z)
            # Atome appartenant a un DOPC
            else:
                dopc[i].append(z)
    return(h2o, dopc, dope)

def make_slices(nb_slices):
    z_range = max_z / nb_slices
    z_list = [min_z]
    z = z_range
    for i in range(nb_slices - 1):
        z += z_range
        z_list.append(z)
    return(z_list)

def number_density(nb_slices, h2o, dopc, dope):
    new_h20 = []
    new_dope = []
    new_dopc = []
    for i in range(nb_slices):
        nb_h20 = 0
        nb_dope = 0
        nb_dopc = 0

        if(i in h2o):
            nb_h20 = len(h2o[i])
        if(i in dope):
            nb_dope = len(dope[i])
        if(i in dopc):
            nb_dopc = len(dopc[i])

        total = nb_h20 + nb_dopc + nb_dope

        if(total == 0):
            new_h20.append(0)
            new_dope.append(0)
            new_dopc.append(0)
        else:
            new_h20.append(round(nb_h20 / total, 4) * 100)
            new_dope.append(round(nb_dope / total, 4) * 100)
            new_dopc.append(round(nb_dopc / total, 4) * 100)
    return(new_h20, new_dopc, new_dope)

def draw_density_profile(h2o, dopc, dope):

    # Gestion si aucuns atomes dans une tranche --> None on ne trace pas de point
    for i in range(len(slices)):
        if(h2o[i] == 0 and dopc[i] == 0 and dope[i] == 0):
            h2o[i] = None
            dopc[i] = None
            dope[i] = None

    y = np.array(slices)
    plt.plot(h2o, y, label="H2O")
    plt.plot(dopc, y, label="DOPC")
    plt.plot(dope, y, label="DOPE")

    plt.xlabel("Densite de nombre (%)")
    plt.ylabel("Distance sur l'axe des z de la membrane (nm)")

    plt.legend()
    plt.show()

frame_id = 5
list_z = traj.xyz[frame_id, : , 2]
min_z = min(list_z)
max_z = max(list_z)

#nb_slices = input("Entrez un nombre pour la division de la membrane en plusieurs tranches et le calcul du profil de densité : ")
slices = make_slices(51)
h2o, dopc, dope = sort_atoms_inside_slices(slices, frame_id)
h2o, dopc, dope = number_density(51, h2o, dopc, dope)

draw_density_profile(h2o, dopc, dope)