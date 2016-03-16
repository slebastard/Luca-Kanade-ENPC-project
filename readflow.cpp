#include "readflow.hpp"
#include "flowcode/Image.h"
#include "flowcode/flowIO.h"



// Seules methode utilisees car les autres contiennent des erreurs. Nous avons donc opte pour utiliser
// 	la fonction de chargement de ground_truth map a partir de fichier .flo fournie sur MiddleBurry


Image<FVector<double, 2>, 2 > flow_from_file(string& path){
	
    // Reads image using flowIO code
	CImageOf<double> img;
	ReadFlowFile(img, path.c_str());
	// Then convert it to Image<FVector<double, 2>, 2 >
	CShape shape = img.Shape();
	Image<FVector<double, 2>, 2 > V(shape.width, shape.height);
	for (int i = 0; i < shape.width; i++){
			for (int j = 0; j < shape.height; j++){
				V(i,j)[0] = (img.Pixel(i,j,0) < 100.0) ? img.Pixel(i,j,0): 0.0; // take out overflow values
				V(i,j)[1] = (img.Pixel(i,j,1) < 100.0) ? img.Pixel(i,j,1): 0.0; // take out overflow values
		}
	}
	return V;
}



Image<double, 2> error_map(const Image<FVector<double, 2>, 2>& ground_truth, const Image<FVector<double, 2>, 2>& flow_estimation, string mode)
{
   
    int w = ground_truth.width();
    int h = ground_truth.height();
    
    if (flow_estimation.width() != w || flow_estimation.height() != h)
        throw ( "Les dimensions des images a comparer ne correspondent pas" );
    
    
    Image<double, 2> error_map;
    if (mode == "HRZT")
        error_map = dim_wise_error(ground_truth, flow_estimation, 0);
    else if (mode == "VERT")
        error_map = dim_wise_error(ground_truth, flow_estimation, 1);
    else if (mode == "NORM")
        error_map = norm_error(ground_truth, flow_estimation);
    else
        throw ( "Methode renseignee inconnue" );
    
    
    return error_map;
}


 Image<double, 2> norm_error(const Image<FVector<double, 2>, 2>& ground_truth, const Image<FVector<double, 2>, 2>& flow_estimation)
 {
 	
 	int w = ground_truth.width();
 	int h = ground_truth.height();

 	if (flow_estimation.width() != w || flow_estimation.height() != h)
 		throw ( "Les dimensions des images a comparer ne correspondent pas" );


 	Image<double, 2> error(w,h);
 	for (int i = 0; i < w; i++)
 		for (int j = 0; j < h; j++)
 			error(i, j) = norm(ground_truth(i, j) - flow_estimation(i, j));

 	return error;
 }



 Image<double, 2> dim_wise_error(const Image<FVector<double, 2>, 2>& ground_truth, const Image<FVector<double, 2>, 2>& flow_estimation, int dim)
 {
 
 	int w = ground_truth.width();
 	int h = ground_truth.height();

 	if (flow_estimation.width() != w || flow_estimation.height() != h)
 		throw ( "Les dimensions des images a comparer ne correspondent pas" );


 	Image<double, 2> error(w,h);
 	if (dim < 0 && dim > 1)
 		throw ( "Dimension de reference invalide pour le calcul de l'erreur" );

 	for (int i = 0; i < w; i++)
 		for (int j = 0; j < h; j++)
 			error(i, j) = abs(ground_truth(i, j)[dim] - flow_estimation(i, j)[dim]);
 	return error;
 }



//  int hex2dec(string hex_bytes)
//  {
//  	/**
//  		hex2dec(hex_bytes):
//  		Transforme un groupement hexadecimal de bytes
//  		$hex_byte$ en decimal
//  	*/

//  	int dec_bytes;
//  	string hex[4];
//  	int dec[4];
//  	std::stringstream stream;
//  	hex[0] = hex_bytes.substr(0, 2);
//  	hex[1] = hex_bytes.substr(2, 2);
//  	hex[2] = hex_bytes.substr(4, 2);
//  	hex[3] = hex_bytes.substr(6, 2);

//  	stream << hex[0];
//  	stream >> std::hex >> dec[0];
//  	stream.str();
//  	stream.clear();
//  	stream << hex[1];
//  	stream >> std::hex >> dec[1];
//  	stream.str();
//  	stream.clear();
//  	stream << hex[2];
//  	stream >> std::hex >> dec[2];
//  	stream.str();
//  	stream.clear();
//  	stream << hex[3];
//  	stream >> std::hex >> dec[3];

//  	dec_bytes = dec[0] * 1 + dec[1] * pow(16, 2) + dec[2] * pow(16, 4) + dec[3] * pow(16, 6);

//  	return dec_bytes;
//  }

// bool is_hexadecimal(string message_bytes)
// {
// 	/**
// 		is_hexadecimal(message_bytes):
// 		*/
// 	// A AFFINER
// 	return true;
// }

// void expand_to_8digits(string& message_byte)
// {
//     //message_byte = std::string(8-message_byte.length(), '0') + message_byte;
//     message_byte = message_byte + std::string(8-message_byte.length(), '0');
// }

// bool is_end_message(string message_byte)
// {
// 	return (message_byte.length() == 2);
// 	// A AFFINER
// }

// double get_bytes_double(istream& stream, int nb_byte)
// {
// 	/**
// 		get_bytes(stream, nb_byte = 4):
// 		Lit les $nb_byte$ prochain bytes du fichier lié
// 		à $stream$ en partant de la position actuelle du
// 		curseur de lecture
// 	*/
// 	int feed_turns = nb_byte / 2;
// 	string current_byte, message_byte = "";

// 	for (int turn = 0; turn < feed_turns; turn++){
// 		stream >> current_byte;
// 		message_byte += current_byte;
// 	}
// 	// SORTIE DE DEBOGUAGE
// 	// A decommenter uniquement pour le deboguage
// 	// cout << "Digit " << stream.tellg() << ": " << message_byte << endl;
// 	if (is_end_message(message_byte))
// 		return -1;

// 	expand_to_8digits(message_byte);
// 	if (!is_hexadecimal(message_byte))
// 		throw ( "Le contenu a traduire ne correspond pas a de l'hexadecimal" );
		

// 	return hex2double(message_byte);
// }


// int get_bytes_int(istream& stream, int nb_byte)
// {
// 	/**
// 		get_bytes(stream, nb_byte = 4):
// 		Lit les $nb_byte$ prochain bytes du fichier lié
// 		à $stream$ en partant de la position actuelle du
// 		curseur de lecture
// 	*/
// 	int feed_turns = nb_byte / 2;
// 	string current_byte, message_byte = "";

// 	for (int turn = 0; turn < feed_turns; turn++){
// 		stream >> current_byte;
// 		message_byte += current_byte;
// 	}
// 	// SORTIE DE DEBOGUAGE
// 	// A decommenter uniquement pour le deboguage
// 	// cout << "Digit " << stream.tellg() << ": " << message_byte << endl;
// 	if (is_end_message(message_byte))
// 		return -1;

// 	expand_to_8digits(message_byte);
// 	if (!is_hexadecimal(message_byte))
// 		throw ( "Le contenu a traduire ne correspond pas a de l'hexadecimal" );
		

// 	return hex2dec(message_byte);
// }

// // hexadecimal to double convertor
// double hex2double(string msg){
//   uint32_t num;
//   sscanf(msg.c_str(), "%x", &num);
//   return *((double*)&num);
// }


// // hexadecimal to double convertor
// uint32_t hex2int(string msg){
//   uint32_t num;
//   sscanf(msg.c_str(), "%x", &num);
//   return num;
// }


// bool right_stream_pos(istream& stream, int row, int col, int width, int height)
// {
// 	/**
// 		right_stream_pos(stream, row, col, width, height):
// 		- Verifie que la position du curseur associe a $stream$
// 		est coherente avec la position dans l'image ($row$
// 		et $col$ renseignees)
// 		- Verifie que la position dans l'image satisfait les
// 		dimensions $width$ et $height$ de l'image
// 	*/
// 	//int stream_pos = stream.tellg();
// 	//bool right_pos = (stream_pos == 30 + 16 * (width*row + col));
// 	//bool inbound = (stream_pos < 12 + 16 * (width*(height - 1) + width));
// 	//return (right_pos && inbound);
// 	return true;
// }
//
//
//
// Image<FVector<double, 2>, 2 > flow_from_txt(string& path)
// {
// 	/**
// 		flow_from_flo(path):
// 		Lit le fichier .flo indique par le chemin $path$,
// 		fabrique et retourne une carte de gradient theorique
// 		a partir du code hexadecimal du fichier. Pour obtenir
// 		les conventions de codage du fichier .flo, se referrer
// 		a la documentation Middlebury
// 	*/

// 	ifstream flow_stream(path);
// 	if (flow_stream)
// 	{
// 		int w, h;
// 		flow_stream.seekg(9, ios::cur);
// 		w = (int)get_bytes_int(flow_stream);
// 		h = (int)get_bytes_int(flow_stream);
// 		Image<FVector<double, 2>, 2 > V(w, h);
// 		//int nb_bytes = (w)* (h)* 8;

// 		for (int row = 0; row < h; row++)
// 		{
// 			for (int col = 0; col < w; col++)
// 			{
// 				if (!right_stream_pos(flow_stream, row, col, w, h))	// On vŽrifie qu'on lit bien les bytes associŽs ˆ ce pixel
// 					throw ( "Erreur indicielle durant la lecture du fichier .flo" );
					
// 				// SORTIE DE DEBOGUAGE
// 				// A decommenter uniquement pour le deboguage
// 				// cout << "Ligne: " << row << " - Colonne: " << col << endl;

// 				cout << get_bytes_int(flow_stream) << " ";
//                 cout << get_bytes_int(flow_stream) << endl;;
//                 //DEBUG
//                 //cout  << ios::cur << " " << V(col, row)[0] << " " << V(col, row)[1] << " ";
// 			}
// 		}
// 		return V;
// 	}
// 	else
// 	{
// 		cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
// 		throw ( "Le chemin renseigne n'a pas ete trouve" );
// 	}
// }

