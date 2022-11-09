
#SCRIPT POUR LE CALCUL DU POTENTIEL ELECTROSTATIQUE TOUT ATOME (ET C ALPHA)


#importation des libraires necessaires
import MDAnalysis as mda
from scipy.spatial import distance
import matplotlib.pyplot as plt
import sys, os, math

#Importation des fichiers et choix de l'utilisateur

def verif():

	'''La fonction verif:
	- Prend comme input: un fichier topologie (format ".tpr"),
	un fichier trajectoire (format ".trr") et figure_name (type str).
	- Permet de vérifier si tous les arguments ont été entrés et répondent aux conditions.
	- Donne comme output: les input s'ils sont corrects et tous spécifiés
	(topol_file, traj_file, figure_name), un message d'erreur explicite sinon.'''

	try:
		topol_file = sys.argv[1]
		traj_file = sys.argv[2]
		figure_name = sys.argv[3]

	#Verification des fichiers:

	#Si biens specifies
	except:
		sys.exit("Erreur: nombre d'options entré incorrect. Le(s) fichier(s) et/ou nom choisi pour le graphe n'ont pas été spécifié(s).\
	\nPour lancer ce programme, il faut entrer comme arguments et dans l'ordre: le fichier topologie, le fichier trajectoire\
 puis le nom désiré pour le fichier en sortie contenant le graphique.")

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

	return topol_file, traj_file, figure_name

#print(verif.__doc__)


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


def calc_frame(topol_file,traj_file):

	'''La fonction calc_frame:
	- Prend comme input: le fichier topologie (format ".tpr")
	et le fichier trajectoire (format ".trr").
	- Permet de parcourir les fichiers frame par frame,
	de trouver et d'envoyer les valeurs des variables correspondantes
	nécaissaires à la fonction potentiel_elect.
	- Donne comme output: le potentiel electrostatique pour chaque frame
	stocké dans potelec_byframe (type list).'''

	#Importation de l'univers MDAnalysis etudie ici
	u = mda.Universe(topol_file, traj_file)

	#Selection des atomes de la proteine et du solvant
	#protein = u.select_atoms("name CA")
	protein = u.select_atoms("protein")
	solvant = u.select_atoms("name OW or name HW1 or name HW2")

	#Selection des charges des atomes de la proteine et du solvant
	prot_charges = list(protein.charges)
	sol_charges = list(solvant.charges)

	#Selection des coordonnées des atomes
	prot_coord = list(protein.positions*0.1)
	sol_coord = list(solvant.positions*0.1)

	#Calcul de chaque distance entre atome proteine et atome solvant
	distance1 = distance.cdist(prot_coord, sol_coord, 'euclidean')

	#Variables pour la boucle
	i=0
	potelec_byframe = []

	for ts in u.trajectory:
		for j in range(len(distance1[0])):
			frame_pot = 0
			r = distance1[i][j]
			#cut-off
			#if r < 3.25:
			if r < 4:
				qi = prot_charges[i]
				qj = sol_charges[i]
				add = potentiel_elect (qi,qj,r)
				frame_pot += add

		potelec_byframe.append(frame_pot)
		i += 1

	return potelec_byframe

#print(calc_frame.__doc__)


def affichage(potelec_byframe, figure_name):

	'''La fonction affichage:
	- Prend comme input: potelec_byframe contenant le potentiel electrostatique
	pour chaque frame (type list) et figure_name (type str).
	- Permet de tracer le graphique des potentiels electrostatiques en fonction
	des frames.
	- Donne comme output: un fichier contenant le graphique représentant
	les potentiels electrostatiques en fonction des frames (format .png)'''

	x = range(0,51)
	y = potelec_byframe
	plt.plot(x,y)
	plt.xlabel('Frames')
	plt.ylabel('Potentiels électrostatiques')
	plt.title("Potentiel électrostatique en fonction des frames")
	plt.savefig(figure_name)

print(affichage.__doc__)


def main():
	topol_file,traj_file,figure_name = verif()
	potelec_byframe = calc_frame(topol_file,traj_file)
	affichage(potelec_byframe,figure_name)

if __name__ == '__main__':
	main()