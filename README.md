Luca-Kanade-ENPC-project
===========================

Méthodes pour le calcul de flot optique - École des Ponts Paristech
L. Cressot & S. Lebastard
Enseignant : P. Monasse


Description du code
===================


Main.cpp
--------
Le main permet de passer des arguments en ligne de commande. Pour connaitres les options, il suffit d'appeller l'aide
en lançant le programme sans options ou bien avec l'option -h.


Une fois les arguments lus sur la ligne de commande, le main charge les images présentent dans un dossier spécifié, ainsi que la map ground_truth si spécifiée, les prétraite si nécéssaire (réduction des dimensions des images), puis effectue les calculs de flot optique.

Si l'option de test est spécifiée, des calculs sont lancés en boucle sur les deux première images pour des paramètres différents et pour les différentes méthodes implémentées, puis les résultats d'erreur (en comparaison avec la ground truth) sont enregistré en format csv.

Sinon, il faut spécifier une méthode de calcul de flot optique et les calculs sont lancés en prenant deux à deux pour toutes les images chargées. Les résultats peuvent être affichés et/ou  enregistrés (voir les options).


Exemples de lancement du programme :

./otpflow -d ../data/images_inputs/ -o ../data/images_outputs/ -m "LK 7" -s -r 50000 
( Charge les images dans ../data/images_inputs/, les traitent 2 à 2 avec la méthode Lucas Kanade avec taille de fenêtre 7, sauvegarde les images dans ../data/images_outputs/ et réduit les images à la résolution 50000 si elles la dépasse)

./optflow -d ../data/images_inputs/ -o ../data/images_outputs/ -t
( Lance une procédure de test sur les deux première images trouvées dans ../data/images_inputs/ et écrit les résultats .csv dans ../data/images_outputs/)



Optflow.cpp
-----------
Ce fichier contient les méthodes de calcul du flot optique ainsi que toutes les fonctions utiles à ces méthodes.



Readflow.cpp
------------
Ce code est destiné à la lecture de fichier .flo en hexadécimal convertis préalablement en .txt, mais il comporte des erreurs et est actuellement non opérationnel.
Pour le remplacer est utilisée une fonction de lecture de fichier .flo (donnée sur MidlleBurry.com) inclue dans le sous-dossier flowcode.
