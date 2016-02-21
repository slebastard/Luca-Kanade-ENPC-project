/*
Projet MOPSI : Flux optique
Loic Cressot  &  Simon Lebastard
Encadrant : Pascal Monasse
janvier 2016
*/

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <string>
#include <vector>
#include <cfloat>

#include <Imagine/Images.h>
#include "Imagine/Common.h"
#include "Imagine/LinAlg.h"

#include "optflow.hpp"

using namespace std;
using namespace Imagine;

// Erreur d'accès dans un Array ici
Image<FVector<float,2> ,2 > flow_Lucas_Kanade(Image<FVector<float,3> >& I1, Image<FVector<float,3> >& I2, int taille_fenetre){

	int w = I1.width(), h = I1.height();
	int w2 = I2.width(), h2 = I2.height();

	if(w!=w2 || h!=h2){
		// Erreur si images de taille différentes
		throw string("Erreur : Images de tailles différentes");
	}

	// 1) Caclcul du gradient
	// On decompose l'image 1 en trois images selon ses 3 composantes R G B 
	Image<float, 2> I1_R(w,h);
	Image<float, 2> I1_G(w,h);
	Image<float, 2> I1_B(w,h);

	for (int i=0; i < w ; i++){
		for (int j=0; j < h ; j++){ 
			I1_R(i,j) = I1(i,j)[0];
			I1_G(i,j) = I1(i,j)[1];
			I1_B(i,j) = I1(i,j)[2];
		}
	}

	// Calcul des composantes des images gradient des 3 composantes R G B
	Image<FVector<float,2>, 2 > gradU_R = image_gradient(I1_R);
	Image<FVector<float,2>, 2 > gradU_G = image_gradient(I1_G);
	Image<FVector<float,2>, 2 > gradU_B = image_gradient(I1_B);

	// 2) Calcul de dtu
	// On calcule maintenant dtu pour les 3 composantes
	Image<float, 2> dtu_R = image_dtu(I1,I2,0);
	Image<float, 2> dtu_G = image_dtu(I1,I2,1);
	Image<float, 2> dtu_B = image_dtu(I1,I2,2);

	//3) Caclcul du flow optique par régression linéaire
	//On crée l'image qui contiendra le flow optique que l'on va retourner
	Image<FVector<float,2> ,2 > V(w,h); // sera retourné une fois rempli

	// Calcul flow optique :
	int n = (int)taille_fenetre/2;
    //const int n_equations = taille_fenetre * taille_fenetre * 3;
	FMatrix<float, 147, 2> A(0.0);
	FVector<float, 147> b(0.0);
	//Matrix<float> A(n_equations,2);
	//Vector<float> b(n_equations);

	// Boucle sur les pixels
	for (int i=0; i < w ; i++){
		for (int j=0; j < h ; j++){
			// On met à 0 la matrice et le vecteur
			A.fill(0.0);
			b.fill(0.0);

			// Boucle sur la fenêtre centrée en (i,j)
			int ind = 0;
			for (int k=max(0,i-n); k <= min(w-1,i+n); k++){
				for (int l=max(0,j-n); l <= min(h-1,j+n); l++){
					//On remplit la matrice et le vecteur
					A(ind,0) = gradU_R(k,l)[0];
					A(ind,1) = gradU_R(k,l)[1];
					b[ind] = -1.0 * dtu_R(k,l);
					ind++;
					A(ind,0) = gradU_G(k,l)[0];
					A(ind,1) = gradU_G(k,l)[1];
					b[ind] = -1.0 * dtu_G(k,l);
					ind++;
					A(ind,0) = gradU_B(k,l)[0];
					A(ind,1) = gradU_B(k,l)[1];
					b[ind] = -1.0 * dtu_B(k,l);
					ind++;
				}
			}
			// On remplit V
			V(i,j) = linSolve(A,b);
		}
	}

	// On retourne le flow optique
	return V;
}

Image<FVector<float, 2>, 2 > flow_Horn_Schunk(Image<FVector<float, 3> >& I1, Image<FVector<float, 3> >& I2, float smoothness, float stop, int iter_max){

	int w = I1.width(), h = I1.height();
	int w2 = I2.width(), h2 = I2.height();

	if (w != w2 || h != h2){
		// Erreur si images de taille différentes
		throw string("Erreur : Images de tailles différentes");
	}

	// 1) Calcul du gradient
	// On decompose l'image 1 en trois images selon ses 3 composantes R G B 
	Image<float, 2> I1_R(w, h);
	Image<float, 2> I1_G(w, h);
	Image<float, 2> I1_B(w, h);

	for (int i = 0; i < w; i++){
		for (int j = 0; j < h; j++){
			I1_R(i, j) = I1(i, j)[0];
			I1_G(i, j) = I1(i, j)[1];
			I1_B(i, j) = I1(i, j)[2];
		}
	}

	// Calcul des composantes des images gradient des 3 composantes R G B
	Image<FVector<float, 2>, 2 > gradU_R = image_gradient(I1_R);
	Image<FVector<float, 2>, 2 > gradU_G = image_gradient(I1_G);
	Image<FVector<float, 2>, 2 > gradU_B = image_gradient(I1_B);

	// 2) Calcul de dtu
	// On calcule maintenant dtu pour les 3 composantes
	Image<float, 2> dtu_R = image_dtu(I1, I2, 0);
	Image<float, 2> dtu_G = image_dtu(I1, I2, 1);
	Image<float, 2> dtu_B = image_dtu(I1, I2, 2);

	/*   Calcul du flow optique par itération sur le modèle de Horn & Schunk    */

	//3.1) Créer une carte de vitesses initiales, idéalement en prenant comme entrée le résultat d'une autre méthode
	//     Pour l'instant l'entrée en vecteurs vitesses sera aléatoire
	Image<FVector<float, 2>, 2 > V_R = init_map(w, h);
	Image<FVector<float, 2>, 2 > V_G = init_map(w, h);
	Image<FVector<float, 2>, 2 > V_B = init_map(w, h);
	Image<FVector<float, 2>, 2 > V_R_ant(w, h);
	Image<FVector<float, 2>, 2 > V_G_ant(w, h);
	Image<FVector<float, 2>, 2 > V_B_ant(w, h);
	Image<FVector<float, 2>, 2 > V_moy_R(w, h);
	Image<FVector<float, 2>, 2 > V_moy_G(w, h);
	Image<FVector<float, 2>, 2 > V_moy_B(w, h);
	bool must_stop = false;
	int iter = 0;
	//3.2) Fixer un critère d'arrêt sur le calcul du gradient
	while (!must_stop && iter < iter_max)
	{
		must_stop = true;
		for (int i = 0; i < w; i++){
			for (int j = 0; j < h; j++){
				V_R_ant(i, j)[0] = V_R(i, j)[0];
				V_R_ant(i, j)[1] = V_R(i, j)[1];
				V_moy_R(i, j)[0] = (1 / 12)*(2 * (V_R(i + 1, j)[0] + V_R(i, j + 1)[0] + V_R(i - 1, j)[0] + V_R(i, j - 1)[0]) + V_R(i + 1, j + 1)[0] + V_R(i + 1, j - 1)[0] + V_R(i - 1, j + 1)[0] + V_R(i - 1, j - 1)[0]);
				V_moy_R(i, j)[1] = (1 / 12)*(2 * (V_R(i + 1, j)[1] + V_R(i, j + 1)[1] + V_R(i - 1, j)[1] + V_R(i, j - 1)[1]) + V_R(i + 1, j + 1)[1] + V_R(i + 1, j - 1)[1] + V_R(i - 1, j + 1)[1] + V_R(i - 1, j - 1)[1]);
				V_R(i, j)[0] = V_moy_R(i, j)[0] - gradU_R(i, j)[0] * (gradU_R(i, j)[0] * V_moy_R(i, j)[0] + gradU_R(i, j)[1] * V_moy_R(i, j)[1] + dtu_R(i, j)) / (smoothness*smoothness + gradU_R(i, j)[0] * gradU_R(i, j)[0] + gradU_R(i, j)[1] * gradU_R(i, j)[1]);
				V_R(i, j)[1] = V_moy_R(i, j)[1] - gradU_R(i, j)[1] * (gradU_R(i, j)[0] * V_moy_R(i, j)[0] + gradU_R(i, j)[1] * V_moy_R(i, j)[1] + dtu_R(i, j)) / (smoothness*smoothness + gradU_R(i, j)[0] * gradU_R(i, j)[0] + gradU_R(i, j)[1] * gradU_R(i, j)[1]);
				if (V_R(i, j)[0] - V_R_ant(i, j)[0] > stop || V_R(i, j)[1] - V_R_ant(i, j)[1] > stop)
					must_stop = false;

				V_G_ant(i, j)[0] = V_G(i, j)[0];
				V_G_ant(i, j)[1] = V_G(i, j)[1];
				V_moy_G(i, j)[0] = (1 / 12)*(2 * (V_G(i + 1, j)[0] + V_G(i, j + 1)[0] + V_G(i - 1, j)[0] + V_G(i, j - 1)[0]) + V_G(i + 1, j + 1)[0] + V_G(i + 1, j - 1)[0] + V_G(i - 1, j + 1)[0] + V_G(i - 1, j - 1)[0]);
				V_moy_G(i, j)[1] = (1 / 12)*(2 * (V_G(i + 1, j)[1] + V_G(i, j + 1)[1] + V_G(i - 1, j)[1] + V_G(i, j - 1)[1]) + V_G(i + 1, j + 1)[1] + V_G(i + 1, j - 1)[1] + V_G(i - 1, j + 1)[1] + V_G(i - 1, j - 1)[1]);
				V_G(i, j)[0] = V_moy_G(i, j)[0] - gradU_G(i, j)[0] * (gradU_G(i, j)[0] * V_moy_G(i, j)[0] + gradU_G(i, j)[1] * V_moy_G(i, j)[1] + dtu_G(i, j)) / (smoothness*smoothness + gradU_G(i, j)[0] * gradU_G(i, j)[0] + gradU_G(i, j)[1] * gradU_G(i, j)[1]);
				V_G(i, j)[1] = V_moy_G(i, j)[1] - gradU_G(i, j)[1] * (gradU_G(i, j)[0] * V_moy_G(i, j)[0] + gradU_G(i, j)[1] * V_moy_G(i, j)[1] + dtu_G(i, j)) / (smoothness*smoothness + gradU_G(i, j)[0] * gradU_G(i, j)[0] + gradU_G(i, j)[1] * gradU_G(i, j)[1]);
				if (V_G(i, j)[0] - V_G_ant(i, j)[0] > stop || V_G(i, j)[1] - V_G_ant(i, j)[1] > stop)
					must_stop = false;

				V_B_ant(i, j)[0] = V_B(i, j)[0];
				V_B_ant(i, j)[1] = V_B(i, j)[1];
				V_moy_B(i, j)[0] = (1 / 12)*(2 * (V_B(i + 1, j)[0] + V_B(i, j + 1)[0] + V_B(i - 1, j)[0] + V_B(i, j - 1)[0]) + V_B(i + 1, j + 1)[0] + V_B(i + 1, j - 1)[0] + V_B(i - 1, j + 1)[0] + V_B(i - 1, j - 1)[0]);
				V_moy_B(i, j)[1] = (1 / 12)*(2 * (V_B(i + 1, j)[1] + V_B(i, j + 1)[1] + V_B(i - 1, j)[1] + V_B(i, j - 1)[1]) + V_B(i + 1, j + 1)[1] + V_B(i + 1, j - 1)[1] + V_B(i - 1, j + 1)[1] + V_B(i - 1, j - 1)[1]);
				V_B(i, j)[0] = V_moy_B(i, j)[0] - gradU_B(i, j)[0] * (gradU_B(i, j)[0] * V_moy_B(i, j)[0] + gradU_B(i, j)[1] * V_moy_B(i, j)[1] + dtu_B(i, j)) / (smoothness*smoothness + gradU_B(i, j)[0] * gradU_B(i, j)[0] + gradU_B(i, j)[1] * gradU_B(i, j)[1]);
				V_B(i, j)[1] = V_moy_B(i, j)[1];
				if (V_B(i, j)[0] - V_B_ant(i, j)[0] > stop || V_B(i, j)[1] - V_B_ant(i, j)[1] > stop)
					must_stop = false;
			}
		}
		iter++;
	}
	//3.3) Calculer le gradient a partir du gradient sur chaque teinte. Pour l'instant on prend la moyenne
	Image<FVector<float, 2>, 2 > V(w, h);
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			V(i, j)[0] = (1 / 3)*(V_R(i, j)[0] + V_G(i, j)[0] + V_B(i, j)[0]);
			V(i, j)[1] = (1 / 3)*(V_R(i, j)[1] + V_G(i, j)[1] + V_B(i, j)[1]);
		}
	}
	// On retourne le flow optique
	return V;
}

Image<FVector<float, 2>, 2 > init_map(int width, int height, float v_min, float v_max)
{
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

// Calcul du gradient spatial d'une image 2D
Image<FVector<float,2> , 2> image_gradient(Image<float, 2>& I){

	Image<FVector<float,2> , 2> gradI(I.width(),I.height());

	for (int i=0; i<gradI.width() ; i++){
		for (int j=0; j<gradI.height() ; j++){  // pour chaque pixel, appel de Image::gradient
			gradI(i,j) = gradient( I, Coords<2>(i,j) );
		}
	}
	return gradI;
}


// Calcul de dtu pour une composante donnée entre deux images
Image<float, 2> image_dtu(Image<FVector<float,3> >& I1, Image<FVector<float,3> >& I2, int c){
	
	Image<float, 2> dtu(I1.width(),I1.height());
	
	for (int i=0; i<dtu.width() ; i++){
		for (int j=0; j<dtu.height() ; j++){  // pour chaque pixel, appel de Image::gradient
			dtu(i,j) = I2(i,j)[c] - I1(i,j)[c];
		}
	}
	return dtu;
}

// Transforme le flow optique en image visualisable
Image<Color, 2> make_flow_visible_hsv(Image<FVector<float,2> ,2 >& I){

	int w = I.width(), h = I.height();
	Image<Color, 2> img(w,h);

	// recherche de la norme maximum
	float maximum = -1.0*FLT_MAX;
	for (int i=0; i<w ; i++){
		for (int j=0; j<h ; j++){
			maximum = max(maximum, I(i,j)*I(i,j));
		}
	}

	// normalisation de la norme des vecteurs entre 0 et 1
	I/=sqrt(maximum);


	// reference vector :
	FVector<float,2> ref_vect(0.0,1.0);
	// Maintenant on peut donner une teinte et valeur à chaque pixel de img avec les valeurs de I
	for (int i=0; i<w ; i++){
		for (int j=0; j<h ; j++){  // pour chaque pixel, appel de Image::gradient
			hsv in;
			float norm = sqrt(I(i,j)*I(i,j)); // Okay
			float angle = angle_oriente(I(i,j),ref_vect)*180.0/M_PI + 180.0; // okay
			in.h = angle ;
			in.v = sqrt(norm) ; // increase a little bit the norm for good printing
			in.s = 1.0;
			// On convertit l'espace HSV en RGB avec une saturation constante à 255*0.7
			rgb out = hsv2rgb(in);
			img(i,j) = Color((int)(255*out.r), (int)(255*out.g), (int)(255*out.b) );
		}
	}
	return img;
}

// Transforme le flow optique en image visualisable gris
Image<Color, 2> make_flow_visible_grey(Image<FVector<float,2> ,2 >& I){

	int w = I.width(), h = I.height();
	Image<Color, 2> img(w,h);

	float maximum = -1.0*FLT_MAX;
	float minimum = FLT_MAX;
	for (int i=0; i<w ; i++){
		for (int j=0; j<h ; j++){  // pour chaque pixel, appel de Image::gradient
			maximum = max(maximum, I(i,j)[0]);
			maximum = max(maximum, I(i,j)[1]);
			minimum = min(minimum, I(i,j)[0]);
			minimum = min(minimum, I(i,j)[1]);
		}
	}

	maximum-= minimum;
	I-= FVector<float,2>(minimum,minimum);
	I/=maximum;

	// Maintenant on peut donner une teinte et valeur à chaque pixel de img avec les valeurs de I
	for (int i=0; i<w ; i++){
		for (int j=0; j<h ; j++){  // pour chaque pixel, appel de Image::gradient
			int level = sqrt(I(i,j)*I(i,j))*255;
			// On convertit l'espace HSV en RGB avec une saturation constante à 255*0.7
			img(i,j) = Color(level, level, level);
		}
	}
	return img;
}


// hsv2rgb
rgb hsv2rgb(hsv in)
{
    double      hh, p, q, t, ff;
    long        i;
    rgb         out;

    if(in.s <= 0.0) {       // < is bogus, just shuts up warnings
        out.r = in.v;
        out.g = in.v;
        out.b = in.v;
        return out;
    }
    hh = in.h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.v * (1.0 - in.s);
    q = in.v * (1.0 - (in.s * ff));
    t = in.v * (1.0 - (in.s * (1.0 - ff)));

    switch(i) {
    case 0:
        out.r = in.v;
        out.g = t;
        out.b = p;
        break;
    case 1:
        out.r = q;
        out.g = in.v;
        out.b = p;
        break;
    case 2:
        out.r = p;
        out.g = in.v;
        out.b = t;
        break;
    case 3:
        out.r = p;
        out.g = q;
        out.b = in.v;
        break;
    case 4:
        out.r = t;
        out.g = p;
        out.b = in.v;
        break;
    case 5:
    default:
        out.r = in.v;
        out.g = p;
        out.b = q;
        break;
    }
    return out;
}


// Calcule un angle orienté entre deux vecteurs
float angle_oriente(FVector<float,2>& v1, FVector<float,2>& v2){
	float L = doubleNorm(v1) * doubleNorm(v2);
	if(L < FLT_EPSILON*5)
		return 0.0;
	else{
		float angle = acos(v1*v2/L);
		float signe = sgn(v1^v2);
		return angle * signe;
	}
}

float sgn(float val){
	if(val < 0.0)
    	return -1.0;
    else
    	return 1.0;
}
