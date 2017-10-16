#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
import mdtraj as md
from mdtraj import element
import numpy as np
import operator as op
import matplotlib.pyplot as plt

# Chargement trajectoires/topologie
traj = md.load('C:\Users\cyprien\Desktop\M2\Modelisation\p3_p4_p5_p8\md_200ns_OK.xtc',
               top='C:\Users\cyprien\Desktop\M2\Modelisation\p3_p4_p5_p8\start.pdb')
topology = traj.topology

def sort_atoms_inside_slices(slices, frame_id):
    """
       Fonction qui permet blablablabla
    """

    h2o = {i: [] for i in range(len(slices) + 1)}
    dopc = {i: [] for i in range(len(slices) + 1)}
    dope = {i: [] for i in range(len(slices) + 1)}

    for residue in topology.residues:
        list_atoms = []
        is_dope = False
        for atom in residue.atoms:
            # Si atome d'H et appartient a un DOP --> DOPE
            if atom.element.symbol == 'H' and not atom.residue.is_water:
                is_dope = True
            # Recuperation pos z atome
            list_atoms.append((traj.xyz[frame_id, atom.index, 2], atom.element.symbol))
        # Parcours des tranches pour placer les atomes.
        # Tranche 1 = z1 --> atomes avec z compris entre z0 (min) et z1.
        for z in list_atoms:
            i = 1
            while z[0] > slices[i] and i < len(slices):
                i += 1
            # Atome appartenant à une molecule H2O
            if atom.residue.is_water:
                h2o[i].append(z[1])
            # Atome appartenant a un DOPE
            elif is_dope:
                dope[i].append(z[1])
            # Atome appartenant a un DOPC
            else:
                dopc[i].append(z[1])
    return(h2o, dopc, dope)

def make_slices(nb_slices):
    """
       Fonction qui permet blablablabla
    """

    z_range = BOX_LENGTH / nb_slices
    z_list = []
    z = 0

    for i in range(nb_slices + 1):
        z_list.append(z)
        z += z_range
    return(z_list)

def number_of_atoms(h2o, dopc, dope):
    """
       Fonction qui permet blablablabla
    """

    new_h2o = []
    new_dope = []
    new_dopc = []

    for i in range(nb_slices + 1):
        new_h2o.append(len(h2o[i]))
        new_dope.append(len(dope[i]))
        new_dopc.append(len(dopc[i]))
    return(new_h2o, new_dopc, new_dope)

def mass_density(h2o, dopc, dope):
    """
       Fonction qui permet blablablabla
    """

    new_h2o = []
    new_dope = []
    new_dopc = []
    for i in range(nb_slices + 1):
        if len(h2o[i]) > 0:
            new_h2o.append(density('H', h2o[i].count('H'))
                           + density('O', h2o[i].count('O')))
        else:
            new_h2o.append(0)

        if (len(dope[i]) > 0):
            new_dope.append(density('C', dope[i].count('C'))
                            + density('P', dope[i].count('P'))
                            + density('O', dope[i].count('O'))
                            + density('N', dope[i].count('N'))
                            + density('H', dope[i].count('H')))
        else:
            new_dope.append(0)

        if (len(dopc[i]) > 0):
            new_dopc.append(density('C', dopc[i].count('C'))
                            + density('P', dopc[i].count('P'))
                            + density('O', dopc[i].count('O'))
                            + density('N', dopc[i].count('N')))
        else:
            new_dopc.append(0)
    return(new_h2o, new_dopc, new_dope)

def density(element, nb_atoms):
    n = nb_atoms / NA # Calcul quantité matière
    # Calcul masse
    if element == 'C':
        m = n * MC
    elif element == 'P':
        m = n * MP
    elif element == 'O':
        m = n * MO
    elif element == 'N':
        m = n * MN
    elif element == 'H':
        m = n * MH
    # Calcul densité
    slice_length = BOX_LENGTH / nb_slices
    vol = (slice_length**slice_length) * 1E-27 # Volume tranche en m^-3
    m = m * 1E-3 # Masse en kg
    d = m / vol # Densité en kg/m^-3
    return(d)

def number_density(h2o, dopc, dope):
    """
       Fonction qui permet blablablabla
    """
    new_h2o = []
    new_dopc = []
    new_dope = []

    for i in range(nb_slices + 1):
        total = h2o[i] + dopc[i] + dope[i]
        if(total == 0):
            new_h2o.append(0)
            new_dopc.append(0)
            new_dope.append(0)
        else:
            # Calcul pourcentage nombre d'atomes en fonction total d'atomes
            new_h2o.append(round(h2o[i] / total, 4) * 100)
            new_dopc.append(round(dopc[i] / total, 4) * 100)
            new_dope.append(round(dope[i] / total, 4) * 100)
    return(new_h2o, new_dopc, new_dope)

def draw_density_profile(h2o, dopc, dope, mass_density):
    """
       Fonction qui permet blablablabla
    """
    # Gestion si aucun atome dans une tranche --> None : on ne trace pas de point
    for i in range(nb_slices + 1):
        if(h2o[i] == 0 and dopc[i] == 0 and dope[i] == 0):
            h2o[i] = None
            dopc[i] = None
            dope[i] = None

    y = np.array(slices)
    plt.plot(h2o, y, label="H2O")
    plt.plot(dopc, y, label="DOPC")
    plt.plot(dope, y, label="DOPE")
    plt.xlabel("Densite (kg/m^-3)") if mass_density else plt.xlabel("Densite en nombre d'atomes (%)")
    plt.ylabel("Boite (nm)")
    plt.legend()
    plt.show()


MO = element.oxygen.mass # Masse molaire Oxygène
MH = element.hydrogen.mass # Masse molaire Hydrogène
MP = element.phosphorus.mass # Masse molaire Phosphore
MN = element.nitrogen.mass # Masse molaire Azote
MC = element.carbon.mass # Masse molaire Carbone
NA = 6.022e23 # Constante d'Avogadro

BOX_LENGTH = 9.0 # 9 nm de longueur en x, y et z


frame_nb = 5

#─────────────────────────────────────── MAIN ─────────────────────────────────────────#

run = True

while(run):
    print("\n", 32 * "─" , "MENU" , 32 * "─")
    print("1. Calcul du profil de densité massique")
    print("2. Calcul du profil de densité en fonction du nombre d'atomes")
    print("3. Quitter")
    print(72 * "─")
    choix = raw_input("Entrez votre choix [1-3] => ")

    if choix == '1':
        nb_slices = input("\nNombre de tranches pour la division de la boîte => ")
        h2o_final = [0] * (nb_slices + 1)
        dopc_final = [0] * (nb_slices + 1)
        dope_final = [0] * (nb_slices + 1)
        slices = make_slices(nb_slices)  # Découpage de la boite en tranches

        print("\nConstruction du graphique en cours ...")
        for i in range(frame_nb + 1):
            h2o, dopc, dope = sort_atoms_inside_slices(slices, i)
            h2o, dopc, dope = mass_density(h2o, dopc, dope)
            h2o_final = map(op.add, h2o_final, h2o)
            dopc_final = map(op.add, dopc_final, dopc)
            dope_final = map(op.add, dope_final, dope)
        draw_density_profile(h2o_final, dopc_final, dope_final, True)

    elif choix == '2':
        nb_slices = input("\nNombre de tranches pour la division de la boîte => ")
        h2o_final = [0] * (nb_slices + 1)
        dopc_final = [0] * (nb_slices + 1)
        dope_final = [0] * (nb_slices + 1)
        slices = make_slices(nb_slices)  # Découpage de la boite en tranches

        print("\nConstruction du graphique en cours ...")
        for i in range(frame_nb + 1):
            h2o, dopc, dope = sort_atoms_inside_slices(slices, i)
            h2o, dopc, dope = number_of_atoms(h2o, dopc, dope)
            h2o_final = map(op.add, h2o_final, h2o)
            dopc_final = map(op.add, dopc_final, dopc)
            dope_final = map(op.add, dope_final, dope)
        h2o, dopc, dope = number_density(h2o_final, dopc_final, dope_final)
        draw_density_profile(h2o, dopc, dope, False)

    elif choix == '3':
        run = False

    else:
        print("\nLa saisie ne correspond à aucun menu, réessayez")