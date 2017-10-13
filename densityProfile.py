#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from mdtraj import element
import operator as op

# Chargement trajectoires/topologie
traj = md.load('C:\Users\cyprien\Desktop\M2\Modelisation\p3_p4_p5_p8\md_200ns_OK.xtc',
               top='C:\Users\cyprien\Desktop\M2\Modelisation\p3_p4_p5_p8\start.pdb')
topology = traj.topology

def sort_atoms_inside_slices(slices, frame_id):
    """
       Fonction qui permet blablablabla
    """
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
            list_pos_z_atoms.append((traj.xyz[frame_id, atom.index, 2], atom.element.mass))
        # Parcours des tranches pour placer les atomes dans les bons intervals
        # Regle: indice 1 (tranche 1) de valeur z1 --> contient tous les atomes
        # de la tranche 1 qui auront un z compris entre z0 (le min) et z1.
        for z in list_pos_z_atoms:
            i = 1
            while z[0] > slices[i] and i < len(slices) - 1:
                i += 1
            # Atome appartenant à une molecule H2O
            if(atom.residue.is_water):
                h2o[i].append(z[1])
            # Atome appartenant a un DOPE
            elif(is_dope):
                dope[i].append(z[1])
            # Atome appartenant a un DOPC
            else:
                dopc[i].append(z[1])
    return(h2o, dopc, dope)

def make_slices():
    """
       Fonction qui permet blablablabla
    """
    z_range = 9 / nb_slices
    z_list = [0]
    z = z_range
    for i in range(nb_slices - 1):
        z += z_range
        z_list.append(z)
    return(z_list)

def number_of_atoms(h2o, dopc, dope):
    """
       Fonction qui permet blablablabla
    """
    new_h2o = []
    new_dope = []
    new_dopc = []
    for i in range(nb_slices):
        new_h2o.append(len(h2o[i]))
        new_dope.append(len(dope[i]))
        new_dopc.append(len(dopc[i]))

    return(new_h2o, new_dope, new_dopc)

def mass_of_atoms(h2o, dopc, dope):
    """
       Fonction qui permet blablablabla
    """
    new_h2o = []
    new_dope = []
    new_dopc = []
    for i in range(nb_slices):
        o_mass_h2o, h_mass_h2o = (0,)*2
        if (len(h2o[i]) > 0):
            nb_o_h2o = h2o[i].count(MO)
            nb_h_h2o = h2o[i].count(MH)
            if(nb_h_h2o > 0):
                h_mass_h2o = (nb_h_h2o / NA) * MH
            if (nb_o_h2o > 0):
                o_mass_h2o = (nb_o_h2o / NA) * MO
            new_h2o.append(o_mass_h2o + h_mass_h2o)
        else:
            new_h2o.append(0)

        o_mass_dope, h_mass_dope, n_mass_dope, c_mass_dope, p_mass_dope = (0,)*5
        if (len(dope[i]) > 0):
            nb_o_dope = dope[i].count(MO)
            nb_h_dope = dope[i].count(MH)
            nb_n_dope = dope[i].count(MN)
            nb_c_dope = dope[i].count(MC)
            nb_p_dope = dope[i].count(MP)
            if(nb_c_dope > 0):
                c_mass_dope = (nb_c_dope / NA) * MC
            if (nb_p_dope > 0):
                p_mass_dope = (nb_p_dope / NA) * MP
            if(nb_o_dope > 0):
                o_mass_dope = (nb_o_dope / NA) * MO
            if (nb_n_dope > 0):
                n_mass_dope = (nb_n_dope / NA) * MN
            if (nb_h_dope > 0):
                h_mass_dope = (nb_h_dope / NA) * MH
            new_dope.append(o_mass_dope + h_mass_dope + n_mass_dope + c_mass_dope + p_mass_dope)
        else:
            new_dope.append(0)

        o_mass_dopc, n_mass_dopc, c_mass_dopc, p_mass_dopc = (0,) * 4
        if (len(dopc[i]) > 0):
            nb_o_dopc = dopc[i].count(MO)
            nb_n_dopc = dopc[i].count(MN)
            nb_c_dopc = dopc[i].count(MC)
            nb_p_dopc = dopc[i].count(MP)
            if(nb_c_dopc > 0):
                c_mass_dopc = (nb_c_dopc / NA) * MC
            if (nb_p_dopc > 0):
                p_mass_dopc = (nb_p_dopc / NA) * MP
            if(nb_o_dopc > 0):
                o_mass_dopc = (nb_o_dopc / NA) * MO
            if (nb_n_dopc > 0):
                n_mass_dopc = (nb_n_dopc / NA) * MN
            new_dopc.append(o_mass_dopc + n_mass_dopc + c_mass_dopc + p_mass_dopc)
        else:
            new_dopc.append(0)
    return(new_h2o, new_dope, new_dopc)

def mass_density(h2o, dopc, dope):
    """
       Fonction qui permet blablablabla
    """
    new_h2o = []
    new_dope = []
    new_dopc = []
    for i in range(nb_slices):
        total = h2o[i] + dopc[i] + dope[i]
        if(total == 0):
            new_h2o.append(0)
            new_dope.append(0)
            new_dopc.append(0)
        else:
            new_h2o.append(round(h2o[i] / total, 4) * 100)
            new_dope.append(round(dope[i] / total, 4) * 100)
            new_dopc.append(round(dopc[i] / total, 4) * 100)
    return(new_h2o, new_dopc, new_dope)

def number_density(h2o, dopc, dope):
    """
       Fonction qui permet blablablabla
    """
    new_h2o = []
    new_dope = []
    new_dopc = []
    for i in range(nb_slices):
        total = h2o[i] + dopc[i] + dope[i]
        if(total == 0):
            new_h2o.append(0)
            new_dope.append(0)
            new_dopc.append(0)
        else:
            new_h2o.append(round(h2o[i] / total, 4) * 100)
            new_dope.append(round(dope[i] / total, 4) * 100)
            new_dopc.append(round(dopc[i] / total, 4) * 100)
    return(new_h2o, new_dopc, new_dope)

def draw_density_profile(h2o, dopc, dope, mass_density):
    """
       Fonction qui permet blablablabla
    """
    # Gestion si aucun atome dans une tranche --> None on ne trace pas de point
    for i in range(nb_slices):
        if(h2o[i] == 0 and dopc[i] == 0 and dope[i] == 0):
            h2o[i] = None
            dopc[i] = None
            dope[i] = None

    y = np.array(slices)
    plt.plot(h2o, y, label="H2O")
    plt.plot(dopc, y, label="DOPC")
    plt.plot(dope, y, label="DOPE")
    plt.xlabel("Densite (km/m^-3)") if mass_density else plt.xlabel("Densite en nombre d'atomes (%)")
    plt.ylabel("Boite (nm)")
    plt.legend()
    plt.show()


MO = element.oxygen.mass # Masse molaire Oxygène
MH = element.hydrogen.mass # Masse molaire Hydrogène
MP = element.phosphorus.mass # Masse molaire Phosphore
MN = element.nitrogen.mass # Masse molaire Azote
MC = element.carbon.mass # Masse molaire Carbone
NA = 6.022e23 # Constante d'Avogadro

# list_z = traj.xyz[i, :, 2]
# min_z = min(list_z)
# max_z = max(list_z)



frame_nb = 2

#─────────────────────────────────────── MAIN ─────────────────────────────────────────#

run = True

while(run):
    print("\n", 31 * "─" , "MENU" , 31 * "─")
    print("1. Calcul du profil de densité massique membranaire")
    print("2. Calcul du profil de densité en fonction du nombre d'atomes")
    print("3. Quitter")
    print(69 * "─")
    choix = input("Entrez votre choix [1-3] => ")
    nb_slices = input("\nEntrez un nombre pour la division de la membrane "
                      "en plusieurs tranches et le calcul du profil de densité : ") + 1

    h2o_final = [0] * nb_slices
    dope_final = [0] * nb_slices
    dopc_final = [0] * nb_slices

    slices = make_slices() # Découpage de la boite en tranches

    if choix == 1:
        print("\nConstruction du graphique en cours ...")
        for i in range(frame_nb + 1):
            h2o, dopc, dope = sort_atoms_inside_slices(slices, i)
            h2o, dopc, dope = mass_of_atoms(h2o, dopc, dope)
            h2o_final = map(op.add, h2o_final, h2o)
            dopc_final = map(op.add, dopc_final, dopc)
            dope_final = map(op.add, dope_final, dope)

        # h2o, dopc, dope = number_density(h2o_final, dopc_final, dope_final)
        draw_density_profile(h2o, dopc, dope, True)
    elif choix == 2:
        print("\nConstruction du graphique en cours ...")
        for i in range(frame_nb + 1):
            h2o, dopc, dope = sort_atoms_inside_slices(slices, i)
            h2o, dopc, dope = number_of_atoms(h2o, dopc, dope)
            h2o_final = map(op.add, h2o_final, h2o)
            dopc_final = map(op.add, dopc_final, dopc)
            dope_final = map(op.add, dope_final, dope)
            h2o, dopc, dope = number_density(h2o_final, dopc_final, dope_final)
        draw_density_profile(h2o, dopc, dope, False)
    elif choix == 3:
        run = False
    else:
        print("\nLa saisie ne correspond à aucun menu, réessayez")