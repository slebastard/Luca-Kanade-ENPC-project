#include "readflow.h"

int hex2dec(string hex_bytes)
{
	/*
	hex2dec(hex_bytes):
	Transforme un groupement hexadecimal de bytes
	$hex_byte$ en decimal
	*/

	int dec_bytes;
	std::stringstream stream;

	stream << hex_bytes;
	stream >> std::hex >> dec_bytes;

	return dec_bytes;
}

int get_bytes(istream& stream, int nb_byte)
{
	/*
	get_bytes(stream, nb_byte = 4):
	Lit les $nb_byte$ prochain bytes du fichier lié
	à $stream$ en partant de la position actuelle du
	curseur de lecture
	*/

	int feed_turns = nb_byte / 2;
	int message;
	string current_byte, message_byte = "";

	for (int turn = 0; turn < feed_turns; turn++)
	{
		stream >> current_byte;
		message_byte += current_byte;
	}

	message = hex2dec(message_byte);
	return message;
}

bool right_stream_pos(istream& stream, int row, int col, int width, int height)
{
	/*
	right_stream_pos(stream, row, col, width, height):
	- Verifie que la position du curseur associe a $stream$
	est coherente avec la position dans l'image ($row$
	et $col$ renseignees)
	- Verifie que la position dans l'image satisfait les
	dimensions $width$ et $height$ de l'image
	*/
	int stream_pos = stream.tellg();
	bool right_pos = (stream_pos == 12 + 16 * (width*row + col));
	bool inbound = (stream_pos < 12 + 16 * (width*(height - 1) + width));
	return (right_pos && inbound);
}

Image<FVector<float, 2>, 2 > init_map(int width, int height, float v_min, float v_max)
{
	/*
	init_map(width, height, v_min, v_max):
	Duplicat de la fonction init_map de optflow.cpp
	A SUPPRIMER LORS DE L'INTEGRATION
	*/
	Image<FVector<float, 2>, 2 > V(width, height);
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			V(i, j)[0] = v_min + (v_max - v_min) * doubleRandom();
			V(i, j)[1] = v_min + (v_max - v_min) * doubleRandom();
		}
	}
	return V;
}

Image<FVector<float, 2>, 2 > flow_from_flo(string& path)
{
	/*
	flow_from_flo(path):
	Lit le fichier .flo indique par le chemin $path$,
	fabrique et retourne une carte de gradient theorique
	a partir du code hexadecimal du fichier. Pour obtenir
	les conventions de codage du fichier .flo, se referrer
	a la documentation Middlebury
	*/
	Image<FVector<float, 2>, 2 > V;
	ifstream flow_stream(path);
	if (flow_stream)
	{
		int w, h;
		flow_stream.seekg(4, ios::cur);
		w = get_bytes(flow_stream);
		h = get_bytes(flow_stream);
		int nb_bytes = (w)* (h)* 8;
		for (int row = 0; row < h; row++)
		{
			for (int col = 0; col < w; col++)
			{
				if (!right_stream_pos(flow_stream, row, col, w, h))	// On vérifie qu'on lit bien les bytes associe a ce pixel
					throw exception("Erreur indicielle durant la lecture du fichier .flo");
				for (int dim = 0; dim < 2; dim++)
				{
					V(row, col)[dim] = get_bytes(flow_stream);
				}
			}
		}
		return V;
	}
	else
	{
		cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
		throw exception("Le chemin renseigne n'a pas ete trouve");
	}
}

Image<float, 2> dim_wise_error(Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation, int dim)
{
	/*
	dim_wise_error(ground_truth, flow_estimation, dim):
	Computes an image of error between the $ground_truth$
	and the $flow_estimation$ based on their local difference
	on the $dim$ dimension
	*/
	int w = ground_truth.width();
	int h = ground_truth.height();

	if (flow_estimation.width() != w || flow_estimation.height() != h)
		throw exception("Les dimensions des images a comparer ne correspondent pas");

	Image<float, 2> error;
	if (dim < 0 && dim > 1)
		throw exception("Dimension de reference invalide pour le calcul de l'erreur");
	for (int row = 0; row < h; row++)
		for (int col = 0; col < w; col++)
			error(row, col) = abs(ground_truth(row, col)[dim] - flow_estimation(row, col)[dim]);
	return error;
}

Image<float, 2> norm_error(Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation, bool quadratic)
{
	/*
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
		throw exception("Les dimensions des images a comparer ne correspondent pas");

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
	/*
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
		throw exception("Les dimensions des images a comparer ne correspondent pas");

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
		throw exception("Methode renseignee inconnue");

	return error_map;
}
