#include "readflow.hpp"

int hex2dec(string hex_bytes)
{
	/**
		hex2dec(hex_bytes):
		Transforme un groupement hexadecimal de bytes
		$hex_byte$ en decimal
	*/

	int dec_bytes;
	string hex[4];
	int dec[4];
	std::stringstream stream;
	hex[0] = hex_bytes.substr(0, 2);
	hex[1] = hex_bytes.substr(2, 2);
	hex[2] = hex_bytes.substr(4, 2);
	hex[3] = hex_bytes.substr(6, 2);

	stream << hex[0];
	stream >> std::hex >> dec[0];
	stream.str();
	stream.clear();
	stream << hex[1];
	stream >> std::hex >> dec[1];
	stream.str();
	stream.clear();
	stream << hex[2];
	stream >> std::hex >> dec[2];
	stream.str();
	stream.clear();
	stream << hex[3];
	stream >> std::hex >> dec[3];

	dec_bytes = dec[0] * 1 + dec[1] * pow(16, 2) + dec[2] * pow(16, 4) + dec[3] * pow(16, 6);

	return dec_bytes;
}

bool is_hexadecimal(string message_bytes)
{
	/**
		is_hexadecimal(message_bytes):
		*/
	// A AFFINER
	return true;
}

void expand_to_8digits(string& message_byte)
{
	while (message_byte.length() < 8)
		message_byte = message_byte + "0";
}

bool is_end_message(string message_byte)
{
	return (message_byte.length() == 2);
	// A AFFINER
}

float get_bytes(istream& stream, int nb_byte)
{
	/**
		get_bytes(stream, nb_byte = 4):
		Lit les $nb_byte$ prochain bytes du fichier li�
		� $stream$ en partant de la position actuelle du
		curseur de lecture
	*/
	int feed_turns = nb_byte / 2;
	string current_byte, message_byte = "";

	for (int turn = 0; turn < feed_turns; turn++){
		stream >> current_byte;
		message_byte += current_byte;
	}
	// SORTIE DE DEBOGUAGE
	// A decommenter uniquement pour le deboguage
	// cout << "Digit " << stream.tellg() << ": " << message_byte << endl;
	if (is_end_message(message_byte))
		return -1;

	expand_to_8digits(message_byte);
	if (!is_hexadecimal(message_byte))
		throw ( "Le contenu a traduire ne correspond pas a de l'hexadecimal" );
		

	return hex2float(message_byte);
}

// hexadecimal to float convertor
float hex2float(string msg){
  uint32_t num;
  sscanf(msg.c_str(), "%x", &num);
  return *((float*)&num);
}

bool right_stream_pos(istream& stream, int row, int col, int width, int height)
{
	/**
		right_stream_pos(stream, row, col, width, height):
		- Verifie que la position du curseur associe a $stream$
		est coherente avec la position dans l'image ($row$
		et $col$ renseignees)
		- Verifie que la position dans l'image satisfait les
		dimensions $width$ et $height$ de l'image
	*/
	//int stream_pos = stream.tellg();
	//bool right_pos = (stream_pos == 30 + 16 * (width*row + col));
	//bool inbound = (stream_pos < 12 + 16 * (width*(height - 1) + width));
	//return (right_pos && inbound);
	return true;
}


Image<FVector<float, 2>, 2 > flow_from_txt(string& path)
{
	/**
		flow_from_flo(path):
		Lit le fichier .flo indique par le chemin $path$,
		fabrique et retourne une carte de gradient theorique
		a partir du code hexadecimal du fichier. Pour obtenir
		les conventions de codage du fichier .flo, se referrer
		a la documentation Middlebury
	*/

	ifstream flow_stream(path);
	if (flow_stream)
	{
		int w, h;
		flow_stream.seekg(9, ios::cur);
		w = (int)get_bytes(flow_stream);
		h = (int)get_bytes(flow_stream);
		Image<FVector<float, 2>, 2 > V(w, h);
		int nb_bytes = (w)* (h)* 8;

		for (int row = 0; row < h; row++)
		{
			for (int col = 0; col < w; col++)
			{
				if (!right_stream_pos(flow_stream, row, col, w, h))	// On v�rifie qu'on lit bien les bytes associe a ce pixel
					throw ( "Erreur indicielle durant la lecture du fichier .flo" );
					
				// SORTIE DE DEBOGUAGE
				// A decommenter uniquement pour le deboguage
				// cout << "Ligne: " << row << " - Colonne: " << col << endl;

				V(col, row)[0] = get_bytes(flow_stream);
				V(col, row)[1] = get_bytes(flow_stream);				
			}
		}
		return V;
	}
	else
	{
		cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
		throw ( "Le chemin renseigne n'a pas ete trouve" );
	}
}

Image<float, 2> dim_wise_error(Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation, int dim)
{
	/**
		dim_wise_error(ground_truth, flow_estimation, dim):
		Computes an image of error between the $ground_truth$
		and the $flow_estimation$ based on their local difference
		on the $dim$ dimension
	*/
	int w = ground_truth.width();
	int h = ground_truth.height();

	if (flow_estimation.width() != w || flow_estimation.height() != h)
		throw ( "Les dimensions des images a comparer ne correspondent pas" );
		

	Image<float, 2> error;
	if (dim < 0 && dim > 1)
		throw ( "Dimension de reference invalide pour le calcul de l'erreur" );
		
	for (int row = 0; row < h; row++)
		for (int col = 0; col < w; col++)
			error(row, col) = abs(ground_truth(row, col)[dim] - flow_estimation(row, col)[dim]);
	return error;
}

Image<float, 2> norm_error(Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation, bool quadratic)
{
	/**
		norm_error(ground_truth, flow_estimation, quadratic = false):
		Computes an image of error between the $ground_truth$
		and the $flow_estimation$ based on their local difference
		regarding the euclidian norm.
		If $quadratic$ is set to true then the euclidian norm
		squared is taken into account.
	*/
	int w = ground_truth.width();
	int h = ground_truth.height();

	if (flow_estimation.width() != w || flow_estimation.height() != h)
		throw ( "Les dimensions des images a comparer ne correspondent pas" );
		

	Image<float, 2> error;
	for (int row = 0; row < h; row++)
		for (int col = 0; col < w; col++)
		{
		error(row, col) = norm2(ground_truth(row, col) - flow_estimation(row, col));
		if (!quadratic)
			error(row, col) = sqrt(error(row, col));
		}
	return error;
}

Image<float, 2> error_map(string feature, Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation)
{
	/**
		error_map(feature, ground_truth, flow_estimation):
		Computes an image of error between the $ground_truth$
		and the $flow_estimation$ based on the specified $feature$
		among	HRZT (horizontal component of gradient vectors)
		VERT (vertical component of gradient vectors)
		NORM (euclidian norm of difference)
		NORM2 (square of euclidian norm of difference)
	*/
	int w = ground_truth.width();
	int h = ground_truth.height();

	if (flow_estimation.width() != w || flow_estimation.height() != h)
		throw ( "Les dimensions des images a comparer ne correspondent pas" );
		

	Image<float, 2> error_map;
	if (feature == "HRZT")
		error_map = dim_wise_error(ground_truth, flow_estimation, 0);
	else if (feature == "VERT")
		error_map = dim_wise_error(ground_truth, flow_estimation, 1);
	else if (feature == "NORM")
		error_map = norm_error(ground_truth, flow_estimation, false);
	else if (feature == "NORM2")
		error_map = norm_error(ground_truth, flow_estimation, true);
	else
		throw ( "Methode renseignee inconnue" );
		

	return error_map;
}