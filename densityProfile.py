#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
import mdtraj as md
from mdtraj import element
import numpy as np
import matplotlib.pyplot as plt

def sort_atoms_inside_slices(slices):
    """
       Fonction qui traite les frames une à une et qui parcours chaque atome de
       chaque molécule pour les ranger dans les différentes tranches.
       @param: slices Tableau contenant les positions z des différentes tranches.
       @return: h2o, dopc, dope Dictionnaires associant les tranches et les symboles
       des atomes pour les 3 molécules étudiées.
    """
    h2o = {i: [] for i in range(len(slices) + 1)}
    dopc = {i: [] for i in range(len(slices) + 1)}
    dope = {i: [] for i in range(len(slices) + 1)}

    # Pour chaque frame
    for f in range(frame_nb + 1):
        print("\nTraitement de la frame : " + str(f))
        # Pour chaque molécule
        for residue in topology.residues:
            list_atoms = []
            is_dope = False
            # Pour chaque atome de la molecule courante
            for atom in residue.atoms:
                # Si atome d'H et appartient a un DOP --> DOPE
                if atom.element.symbol == 'H' and not atom.residue.is_water:
                    is_dope = True
                # Recuperation pos z atome et symbole de l'élément
                list_atoms.append((traj.xyz[f, atom.index, 2], atom.element.symbol))
            # Parcours des tranches pour placer les atomes.
            # Tranche 1 : z1 --> contient les atomes avec z compris entre z0 (min) et z1.
            for z in list_atoms:
                i = 1
                while z[0] > slices[i] and i < len(slices)- 1:
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
       Fonction qui permet de diviser la longueur de la boîte par le
       nombre de tranches souhaitées.
       @param: nb_slices Le nombre de tranches défini par l'utilisateur
       @return: z_list Les positions z des tranches
    """
    z_range = Z_BOX / nb_slices
    z_list = []
    z = 0
    for i in range(nb_slices + 1):
        z_list.append(z)
        z += z_range
    return(z_list)

def number_of_atoms(h2o, dopc, dope):
    """
       Fonction qui permet de compter le nombre d'atomes de h2o, dopc et dope
       pour chaque tranche.
       @param: h2o, dopc, dope Dictionnaires associant les tranches et les
       symboles des atomes pour les 3 molécules étudiées.
       @return: new_h2o, new_dopc, new_dope Les listes contenant les nombres
       d'atomes pour chaque tranche.
    """
    new_h2o = []
    new_dope = []
    new_dopc = []

    for i in range(nb_slices + 1):
        new_h2o.append(len(h2o[i]))
        new_dope.append(len(dope[i]))
        new_dopc.append(len(dopc[i]))
    return(new_h2o, new_dopc, new_dope)

def mass_of_atoms(h2o, dopc, dope):
    """
       Fonction qui permet de calculer la masse totale des atomes de h2o, dopc, dope
       contenus dans chaque tranche.
       @param: h2o, dopc, dope Dictionnaires associant les tranches et les symboles
       des atomes pour les 3 molécules étudiées.
       @return: new_h2o, new_dopc, new_dope Les listes contenant la masse totale des
       atomes pour chaque molécule pour chaque tranche.
    """
    new_h2o = []
    new_dope = []
    new_dopc = []

    for i in range(nb_slices + 1):
        if len(h2o[i]) > 0:
            new_h2o.append(mass('H', h2o[i].count('H'))
                           + mass('O', h2o[i].count('O')))
        else:
            new_h2o.append(0)

        if (len(dope[i]) > 0):
            new_dope.append(mass('C', dope[i].count('C'))
                            + mass('P', dope[i].count('P'))
                            + mass('O', dope[i].count('O'))
                            + mass('N', dope[i].count('N'))
                            + mass('H', dope[i].count('H')))
        else:
            new_dope.append(0)

        if (len(dopc[i]) > 0):
            new_dopc.append(mass('C', dopc[i].count('C'))
                            + mass('P', dopc[i].count('P'))
                            + mass('O', dopc[i].count('O'))
                            + mass('N', dopc[i].count('N')))
        else:
            new_dopc.append(0)
    return(new_h2o, new_dopc, new_dope)

def mass(element, nb_atoms):
    """
       Fonction qui permet de calculer la masse pour un type d'atome.
       @param: element Le symbole de l'atome
       @param: nb_atoms Le nombre d'atomes pour l'élément concerné
       @return: m La masse ramenée en kg
    """
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
    return(m *1E-3)

def density(h2o, dopc, dope):
    """
       Fonction qui permet de calculer la densité des h2o, dopc et dope
       pour chaque tranches.
       @param: h2o, dopc, dope Les listes contenant les masses totales pour
       chaque tranche.
       @return: h2o, dopc, dope Les listes contenant les densités des
       molécules pour chaque tranche.
    """
    # Calcul volume tranche
    z_slice = Z_BOX / nb_slices
    x_slice = X_BOX / nb_slices
    y_slice = Y_BOX / nb_slices
    vol = (z_slice * x_slice * y_slice) * 1E-27 # Volume tranche en m^-3

    for i in range(nb_slices + 1):
        # Calcul densité en kg/m^-3
        h2o[i] = h2o[i] / vol
        dopc[i] = dopc[i] / vol
        dope[i] = dope[i] / vol
    return(h2o, dopc, dope)

def density_mean(h2o, dopc, dope):
    """
       Fonction qui permet calculer la moyenne des densités en fonction
       du nobmre de frame parcourues.
       @param: h2o, dopc, dope Les listes contenant les densités des
       molécules pour chaque tranche.
       @return: h2o, dopc, dope Les listes contenant les densités des
       molécules moyennées sur le nombre de frames.
    """
    for i in range(nb_slices + 1):
        # Calcul densité en kg/m^-3
        h2o[i] = h2o[i] / (frame_nb + 1)
        dopc[i] = dopc[i] / (frame_nb + 1)
        dope[i] = dope[i] / (frame_nb + 1)
    return (h2o, dopc, dope)

def number_density(h2o, dopc, dope):
    """
       Fonction qui permet calculer la densité en fonction des nombres
       d'atomes des molécules pour chaque tranche.
       @param: h2o, dopc, dope Les listes contenant les nombres d'atomes
       des molécules pour chaque tranche.
       @return: h2o, dopc, dope Les listes contenant les rapports pour chaque
       tranche du nombre d'atomes d'une molecule / nombre total d'atomes des 3
       molécules.
    """
    for i in range(nb_slices + 1):
        total = h2o[i] + dopc[i] + dope[i]
        if(total == 0):
            h2o[i] = 0
            dopc[i] = 0
            dope[i] = 0
        else:
            # Calcul pourcentage nombre d'atomes en fonction total d'atomes
            h2o[i] = round(h2o[i] / total, 4) * 100
            dopc[i] = round(dopc[i] / total, 4) * 100
            dope[i] = round(dope[i] / total, 4) * 100
    return(h2o, dopc, dope)

def draw_density_profile(h2o, dopc, dope, mass_density):
    """
       Méthode qui permet de tracer les profils de densités.
       @param: h2o, dopc, dope Les listes contenant les densités massiques/
       nombres d'atomes des molécules pour chaque tranche.
       @param: mass_density Permet de tracer la courbe en fonction des densités si vrai,
       sinon trace en fonction du nombre d'atomes.
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
    plt.xlabel("Densite (kg/m^-3)") if mass_density \
        else plt.xlabel("Densite en nombre d'atomes (%)")
    plt.ylabel("Boite (nm)")

    plt.legend()
    plt.show()


MO = element.oxygen.mass # Masse molaire Oxygène
MH = element.hydrogen.mass # Masse molaire Hydrogène
MP = element.phosphorus.mass # Masse molaire Phosphore
MN = element.nitrogen.mass # Masse molaire Azote
MC = element.carbon.mass # Masse molaire Carbone
NA = 6.022E23 # Constante d'Avogadro

X_BOX = 9.2042 # 9 nm de longueur en x, y et z
Y_BOX = 9.2042
Z_BOX = 7.8801

# Récupérer toutes les frames : topology.n_frames
# Ici le nombre est défini en brut pour des questions de performances.
frame_nb = 100

# Chargement trajectoires/topologie
xtc = raw_input("Entrez le chemin de votre fichier xtc => ")
pdb = raw_input("\nEntrez le chemin de votre fichier pdb => ")
print("\nChargement des fichiers en cours ...")
traj = md.load(xtc, top=pdb)
topology = traj.topology

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
        slices = make_slices(nb_slices)
        print("\nConstruction du graphique en cours ...")

        h2o, dopc, dope = sort_atoms_inside_slices(slices)
        h2o, dopc, dope = mass_of_atoms(h2o, dopc, dope)
        h2o, dopc, dope = density(h2o, dopc, dope)
        h2o, dopc, dope = density_mean(h2o, dopc, dope)
        draw_density_profile(h2o, dopc, dope, True)

    elif choix == '2':
        nb_slices = input("\nNombre de tranches pour la division de la boîte => ")
        slices = make_slices(nb_slices)
        print("\nConstruction du graphique en cours ...")

        h2o, dopc, dope = sort_atoms_inside_slices(slices)
        h2o, dopc, dope = number_of_atoms(h2o, dopc, dope)
        h2o, dopc, dope = number_density(h2o, dopc, dope)
        draw_density_profile(h2o, dopc, dope, False)

    elif choix == '3':
        run = False

    else:
        print("\nLa saisie ne correspond à aucun menu, réessayez")
