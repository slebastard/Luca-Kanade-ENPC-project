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
* Lit le fichier .flo indique par le chemin $path$, en utilisant une fonction annexe de lecture de fichiers .flo donnee sur Middleburry
*/
Image<FVector<double, 2>, 2 > flow_from_file(string& path);


/**
 * Calcule une map d'erreur entre un flow optique et une ground truth
 */
Image<double, 2> error_map(const Image<FVector<double, 2>, 2>& ground_truth, const Image<FVector<double, 2>, 2>& flow_estimation, string mode = "NORM");

/**
 *  Calcul local d'erreur verticale et horizontale
 */
Image<double, 2> dim_wise_error(const Image<FVector<double, 2>, 2>& ground_truth,const Image<FVector<double, 2>, 2>& flow_estimation, int dim);


/**
 * Calcul local d'erreur L2
 */
Image<double, 2> norm_error(const Image<FVector<double, 2>, 2>& ground_truth, const Image<FVector<double, 2>, 2>& flow_estimation);



// TOUTES LES FONCTIONS QUI SUIVENT NE SONT PAS UTILISEES CAR ELLES CONTIENNENT DES ERREURS
// ========================================================================================
// ========================================================================================
// ========================================================================================
// ========================================================================================
// ========================================================================================
// ========================================================================================


/**
* Transforme un groupement hexadecimal de bytes en int
*/
int hex2dec(string hex_byte);

/**
* Verifie si la string $message_bytes$ correspond a un code hexadecimal valide
*/
bool is_hexadecimal(string message_bytes);

/**
* Etend la string $message_byte$ dans le cas ou celle-ci fait moins de 8 caracteres
*/
void expand_to_8digits(string& message_byte);

/**
* Verifie si le message message_byte recupere correspond a un code de fin de lecture
*/
bool is_end_message(string message_byte);

/**
* Lit les $nb_byte$ prochain bytes du fichier lié à $stream$ en partant de la position actuelle du
* curseur de lecture
*/
double get_bytes_double(istream& stream, int nb_byte = 4);
int get_bytes_int(istream& stream, int nb_byte = 4);

/**
* Convertir un hexadecimal en double
*/
double hex2double(string msg);
uint32_t hex2int(string msg);


/**
*	- Verifie que la position du curseur associe a $stream$ est coherente avec la position dans l'image ($row$ et $col$ renseignees)
*	- Verifie que la position dans l'image satisfait les dimensions $width$ et $height$ de l'image
*/
bool right_stream_pos(istream& stream, int row, int col, int width, int height);

/**
* Lit le fichier .flo indique par le chemin $path$, fabrique et retourne une carte de gradient theorique a partir
* du code hexadecimal du fichier. Pour obtenir les conventions de codage du fichier .flo, se referrer a la documentation Middlebury
*/
Image<FVector<double, 2>, 2 > flow_from_txt(string& path);



#endif