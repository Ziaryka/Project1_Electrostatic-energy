# Project1_Electrostatic-energy

MD_Project9.py est un programme dédié aux calculs d’intéractions électrostatiques entre une protéine et la couche de solvatation lors d’une simulation de dynamique moléculaire.

Ce script se compose en deux parties : la première consiste à l'exécuter suivant un calcul tout atome puis un calcul avec les C_alpha. Quant à la deuxième partie (MD_Project9_resid), celle ci s’applique au calcul des 

potentiels électrostatiques entre les résidus de la protéine et la couche de solvatation. 



##### Installation :  (commune au deux parties) ######

Avant d’exécuter MD_Project9.py, il est impératif de suivre les lignes de commande présentent ci dessous : 

# Activation de l’environnement Conda : (lignes séparées)

conda config --add channels conda-forge

conda update --yes conda

# Installation des packages nécessaires :  MDAnalysis, Numpy, Scipy, Matplotlib.pyplot, Scipy, os et math.

scipy & numpy : sudo apt-get install python-numpy python-scipy

MDAnalysis :

# création de l’environnement : 
conda create --yes -n mdaenv python=3.6 

# installation :
conda install --yes -n mdaenv MDAnalysis MDAnalysisTests 

# activation : 

source activate mdaenv



#### Les fichiers d’entrée (input) : ####

A partir d’un fichier au format .tpr (topol.tpr) utilisé comme topologie de dynamique moléculaire et appelé dans ce script “topol_file” et d’un second fichier au format .trr (md_1BTA_100ns_dt2ns.trr) des trajectoires 

résultantes de la simulation de la protéine Barstar appelé dans ce programme traj_file.

argument 1 = topol.tpr

argument 2 = md_1BTA_100ns_dt2ns.trr

argument 3 = figure_name (choisissez le nom de votre figure exemple : figure_potentiel)

argument 4 =  (utilisé pour la deuxième partie avec les résidus)



#### Vous trouverez dans ce script : ####


Une fonction “verif()”  qui permet de vérifier si tous les arguments ont été entrés et répondent aux conditions. Le but de cette fonction est de s’assurer que les fichiers d’entrée ont bien été importés, ne sont pas vide et présentent la bonne extension.


Une fonction “potentiel_elect()” qui permet de calculer les énergies d’interactions. Elle prend en entrée les charges partielles de la protéine et du solvant (qi,qj) et la distance entre les atomes. Cette fonction renvoie le potentiel électrostatique.


Une fonction “calc_frame()” qui permet de calculer les énergies d’interactions tout au long des frames des trajectoires. Elle prend en input la topologie et la trajectoire (topol_file, traj_file). Cette fonction permet de parcourir les fichiers frames (modèles) de trouver et d’envoyer les valeurs des variables correspondantes nécessaires à la fonction potentiel_elect. Cette fonction renvoie comme fichier de sortie (output) potelec_byframe, une liste contenant le potentiel de chaque frame.


Une fonction affichage() qui permet la visualisation des variations du potentiel électrostatique obtenu lors d’une dynamique moléculaire. Elle prend en input : potelec_byframe, et figure_name. Elle permet de tracer le graphique des potentiels électrostatiques en fonction des frames. Cette fonction renvoie comme output un fichier.png représentant les potentiels électrostatiques.


Une fonction main() qui permet de démarrer et d’exécuter le programme en tant que programme principal. Elle renvoie les sorties de chaque fonction définie précédemment.



#### Exécution du script : ####

Mettre dans un seul répertoire les deux fichiers qui serviront de fichier d'entrée (input), “md_1BTA_100ns_dt2ns.trr” pour les trajectoires et “topol.tpr” pour la topologie de la protéine Barstar.

# Dans un terminal : 

Activer l’environnement Conda.

Lancer python3 

Installer et importer les packages MDAnalysis, Numpy, Scipy, Matplotlib.pyplot, os et math.

Lancer le programme. 

## Ci dessous un exemple de commande pour lancer le programme : ## 

python3 MD_Project9.py topol.tpr md_1BTA_100ns_dt2ns.trr figure_potentiel


#### Recommandations : ####

Afin d’appliquer le cutoff, il vous est recommandé de décommenter la ligne 108 du script.

Décommenter également la ligne 84. Cette ligne vous permettra de lancer le calcul du potentiel électrostatique entre les atomes C_alpha de la protéine et les atomes du solvant. Par conséquent, il vous est suggéré de commenter la sélection des atomes de la protéine (ligne 85).

Décommenter les  __doc__ présents dans chaque fonction afin d’avoir la description de ces dernières dans le terminal lors de l’exécution du script (uniquement la ligne 53 pour la deuxième partie).
