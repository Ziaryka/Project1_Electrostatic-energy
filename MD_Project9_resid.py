#SCRIPT POUR LE CALCUL DU POTENTIEL ELECTROSTATIQUE DES RESIDUS

#importation des libraires necessaires
import MDAnalysis as mda
from scipy.spatial import distance
import matplotlib.pyplot as plt
import numpy as np
import sys, os, math

#Importation des fichiers et choix de l'utilisateur

try:
	topol_file = sys.argv[1]
	traj_file = sys.argv[2]
	figure_name_1 = sys.argv[3]
	figure_name_2 = sys.argv[4]

#Verification des fichiers:

#Si biens specifies
except:
	sys.exit("Erreur: nombre d'options entré incorrect. Le(s) fichier(s) et/ou nom choisi pour le graphe n'ont pas été spécifié(s).\
\nPour lancer ce programme, il faut entrer comme arguments et dans l'ordre: le fichier topologie, le fichier trajectoire\
 puis les noms désirés pour les 2 fichiers en sortie contenant les graphiques.")

#Si bons formats
if not topol_file.endswith(".tpr"):
	sys.exit("Erreur: Mauvais format. Le fichier topologie doit être au format \"tpr\".")
if not traj_file.endswith(".trr"):
	sys.exit("Erreur: Mauvais format. Le fichier trajectoire doit être au format \"trr\".")

#Si existent et bien remplis
if not os.path.exists(topol_file) or os.stat(topol_file).st_size == 0:
	sys.exit("Le fichier topologie spécifié est vide ou il n'existe pas.")
if not os.path.exists(traj_file) or os.stat(traj_file).st_size == 0:
	sys.exit("Le fichier trajectoire spécifié est vide ou il n'existe pas.")


def potentiel_elect (qi,qj,r):

	'''La fonction potentiel_elect:
	- Prend comme input: la charge qi (type float), la charge qj (type float),
	la distance r entre ces deux atomes chargés (type float).
	- Permet de calculer le potentiel electrostatique pour chaque paire d'atomes protéine-solvant.
	- Donne comme output: le potentiel electrostatique "energy" (type float).'''

	pi = math.pi
	epsilon = 80.2

	energy = ((qi*qj)/(4*pi*epsilon*r))
	return energy

#print(potentiel_elect.__doc__)


#Importation de l'univers MDAnalysis etudie ici
u = mda.Universe(topol_file, traj_file)

#Selection des residus de la proteine et des atomes du solvant
prot_res = u.select_atoms("protein and resid 0:88")
protein = u.select_atoms("protein")
solvant = u.select_atoms("name OW or name HW1 or name HW2")

#Selection des charges des residues et du solvant
resid_charges = list(prot_res.residues.charges)
sol_charges = list(solvant.charges)

#Selection des coordonnees des atomes du solvant
sol_coord = list(solvant.positions*0.1)

#Selection des coordonnees du centre de masse des residus de la proteine
resid_com =  np.array([r.atoms.center_of_mass() for r in u.select_atoms("protein").residues])*0.1

#Calcul de chaque distance entre residu proteine et atome solvant
dist = distance.cdist(resid_com, sol_coord, 'euclidean')

#Listes stocakges pour la boucle
pot_resid_frame1 = []
pot_resid_frame2 = []

#Boucle permettant de parcourir les fichiers frame par frame,
#de trouver et d'envoyer les valeurs des variables correspondantes
#nécessaires à la fonction potentiel_elect pour la première et dernière frame.

for ts in u.trajectory: 

	#Pour la premiere frame
	if ts.frame == 0:
		for res in range(len(resid_charges)):
			for sol in range (len(sol_charges)):
				pot = 0
				qi = resid_charges[res]
				qj = sol_charges[sol]
				r = dist[res][sol]
				#cut-off
				if r < 3.25:
					pot += potentiel_elect(qi,qj,r)
			pot_resid_frame1.append(pot)

	#Pour la derniere frame
	if ts.frame == 50:
		for res in range(len(resid_charges)):
			for sol in range (len(sol_charges)):
				pot = 0
				qi = resid_charges[res]
				qj = sol_charges[sol]
				r = dist[res][sol]
				#cut-off
				if r < 4:
					pot += potentiel_elect(qi,qj,r)
			pot_resid_frame2.append(pot)


#Plot des graphiques des potentiels electrostatiques en fonction des frames.
#Donne comme output: un fichier contenant le graphique représentant
#les potentiels electrostatiques en fonction des frames (format .png).

#Plot pour premiere frame
x = range(len(prot_res.residues))
y = pot_resid_frame1
plt.plot(x,y)
plt.xlabel('Résidues')
plt.ylabel('Potentiels électrostatiques')
plt.title("Potentiel électrostatique en fonction des résidues pour la premiere frame")
plt.savefig(figure_name_1)

#Plot pour denriere frame
x = range(len(prot_res.residues))
y = pot_resid_frame2
plt.plot(x,y)
plt.xlabel('Résidues')
plt.ylabel('Potentiels électrostatiques')
plt.title("Potentiel électrostatique en fonction des résidues pour la dernière frame")
plt.savefig(figure_name_2)