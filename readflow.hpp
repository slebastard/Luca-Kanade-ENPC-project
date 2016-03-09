#ifndef READFLOW_H
#define READFLOW_H


#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

#include <Imagine/Images.h>
#include "Imagine/Common.h"
using namespace Imagine;

/**
* Transforme un groupement hexadecimal de bytes
* $hex_byte$ en decimal
* @param hex_byte la chaine de caractere representant lhexadecimal a traduire
* @return nombre decimal traduit
*/
int hex2dec(string hex_byte);

/**
* Verifie si la string $message_bytes$ correspond a un code hexadecimal valide
* @param message_bytes la chaine de caractere representant l'hexadecimal a inspecter
* @return true si le nombre examine correspond a priori a un code hexadecimal valide, false sinon
*/
bool is_hexadecimal(string message_bytes);

/**
* Etend la string $message_byte$ dans le cas ou celle-ci fait moins de 8 caracteres
* @param message_byte la chaine de caractere representant le message hexadecimal a completer par des 0
*/
void expand_to_8digits(string& message_byte);

/**
* Verifie si le message message_byte recupere correspond a un code de fin de lecture
* @param message_byte la chaine de caractere representant l'hexadecimal a examiner
*/
bool is_end_message(string message_byte);

/**
* Lit les $nb_byte$ prochain bytes du fichier lié
* à $stream$ en partant de la position actuelle du
* curseur de lecture
* @param stream doit etre un ifstream pointant sur le fichier sur lequel on va lire les bytes
* @param nb_byte le nombre de bytes a lire durant l'operation. 1 byte hexadecimal = 2 bits hexadecimaux = 256 possibilites
* @return le contenu decimal des bytes lus, si les bytes lus correspondent a de l'hexadecimal
*/
float get_bytes(istream& stream, int nb_byte = 4);

/**
* Convertir un hexadecimal en float
*/
float hex2float(string msg);

/**
* 	right_stream_pos(stream, row, col, width, height):
*	- Verifie que la position du curseur associe a $stream$
*	est coherente avec la position dans l'image ($row$
*	et $col$ renseignees)
*	- Verifie que la position dans l'image satisfait les
*	dimensions $width$ et $height$ de l'image
* @param stream doit etre un ifstream pointant sur le fichier sur lequel on va lire les bytes
* @param row la ligne correspondant au pixel dont on souhaite verifier la correspondance
* @param col la colonne correspondant au pixel dont on souhaite verifier la correspondance
* @param width la largeur de l'image a examiner
* @param height la hauteur de l'image a examiner
* @return true si la position de lecture correspond, false sinon
*/
bool right_stream_pos(istream& stream, int row, int col, int width, int height);

/**
* 	flow_from_flo(path):
* Lit le fichier .flo indique par le chemin $path$,
* fabrique et retourne une carte de gradient theorique
* a partir du code hexadecimal du fichier. Pour obtenir
* les conventions de codage du fichier .flo, se referrer
* a la documentation Middlebury
* @param path est un chemin absolu ou relatif vers le fichier flo a partir duquel on va calculer le ground_truth du gradient
* @return le ground_truth du gradient du flot optique
*/
Image<FVector<float, 2>, 2 > flow_from_txt(string& path);

/**
* 	dim_wise_error(ground_truth, flow_estimation, dim):
* Computes an image of error between the $ground_truth$
* and the $flow_estimation$ based on their local difference
* on the $dim$ dimension
* @param ground_truth est la carte de vitesse reelle selon Middlebury
* @param flow_estimation est l'estimation de l'erreur produite a exteriori
* @param dim est la diemnsions selon laquelle on calcule l'erreur d'estimation
* @see error_map pour l'implementation de cette fonction a la generation d'une carte d'erreur
* @return une carte d'erreur d'estimation
*/
Image<float, 2> dim_wise_error(Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation, int dim);

/**
* 	norm_error(ground_truth, flow_estimation, quadratic):
* Computes an image of error between the $ground_truth$
* and the $flow_estimation$ based on their local difference
* regarding the euclidian norm.
* If $quadratic$ is set to true then the euclidian norm
* squared is taken into account.
* @param ground_truth est la carte de vitesse reelle selon Middlebury
* @param flow_estimation est l'estimation de l'erreur produite a exteriori
* @param quadratic controle la mise au carre de l'erreur. Vaut false par defaut, c'est-a-dire que l'erreur ne sera pas mise au carre
* @see error_map pour l'implementation de cette fonction a la generation d'une carte d'erreur
* @return une carte d'erreur d'estimation
*/
Image<float, 2> norm_error(Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation, bool quadratic = false);

/**
* error_map(feature, ground_truth, flow_estimation):
* Computes an image of error between the $ground_truth$
* and the $flow_estimation$ based on the specified $feature$
* among	HRZT (horizontal component of gradient vectors)
*		VERT (vertical component of gradient vectors)
*		NORM (euclidian norm of difference)
*		NORM2 (square of euclidian norm of difference)
* @feature est le feature selon lequel on va calculer la carte d'erreur:
*		HRZT (horizontal component of gradient vectors)
*		VERT (vertical component of gradient vectors)
*		NORM (euclidian norm of difference)
*		NORM2 (square of euclidian norm of difference)
* @param ground_truth est la carte de vitesse reelle selon Middlebury
* @param flow_estimation est l'estimation de
* @return une carte d'erreur d'estimation
*/
Image<float, 2> error_map(string feature, Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation);



#endif