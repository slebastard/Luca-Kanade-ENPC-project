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
#include <cmath>

#include <Imagine/Images.h>
#include "Imagine/Common.h"
#include "Imagine/LinAlg.h"

#include "optflow.hpp"

using namespace std;
using namespace Imagine;

// Calcul du flow optique par la methode Lucas et Kanade
Image<FVector<float,2> ,2 > flow_Lucas_Kanade(Image<FVector<float,3> >& I1, Image<FVector<float,3> >& I2, int taille_fenetre){

	int w = I1.width(), h = I1.height();
	int w2 = I2.width(), h2 = I2.height();

	if(w!=w2 || h!=h2){
		// Erreur si images de taille différentes
		throw string("Erreur : Images de tailles différentes");
	}

	// 1) Calcul du gradient
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



// Calcul du flow optique par la methode Horn et Schunk itérative
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
		for (int i = 1; i < w-1; i++){
			for (int j = 1; j < h-1; j++){
				// sauvegarder la valeur précédente de V_R
				V_R_ant(i, j) = V_R(i, j);
				// calcul de la moyenne locale R
				V_moy_R(i, j) = V_moy(V_R,i,j);
				// calcul du nouveau V_R
				V_R(i,j) = nouveau_V_c( V_moy_R(i, j), gradU_R(i, j), dtu_R(i,j), smoothness);
				
				if ( maxNorm(V_R(i, j) - V_R_ant(i, j)) > stop )
					must_stop = false;

				// sauvegarder la valeur précédente de V_G
				V_G_ant(i, j) = V_G(i, j);
				// calcul de la moyenne locale G
				V_moy_G(i, j) = V_moy(V_G,i,j);
				// calcul du nouveau V_G
				V_G(i,j) = nouveau_V_c( V_moy_G(i, j), gradU_G(i, j), dtu_G(i,j), smoothness);

				if ( maxNorm(V_G(i, j) - V_G_ant(i, j) ) > stop )
					must_stop = false;

				// sauvegarder la valeur précédente de V_B
				V_B_ant(i, j) = V_B(i, j);
				// calcul de la moyenne locale B
				V_moy_B(i, j) = V_moy(V_B,i,j);
				// calcul du nouveau V_B
				V_B(i,j) = nouveau_V_c( V_moy_B(i, j), gradU_B(i, j), dtu_B(i,j), smoothness);

				if ( maxNorm(V_B(i, j) - V_B_ant(i, j) ) > stop )
					must_stop = false;
			}
		}
		iter++;
	}
    cout << iter-1 << " itérations" << endl;

	//3.3) Calculer le gradient a partir du gradient sur chaque teinte. Pour l'instant on prend la moyenne
	Image<FVector<float, 2>, 2 > V = ( V_R + V_G + V_B ) / 3.0 ;
	
	// On retourne le flow optique
	return V;
}

// Calcul du flow optique par la methode Horn et Schunk Huber L1
Image<FVector<float, 2>, 2 > flow_Horn_Schunk_HuberL1(Image<FVector<float, 3> >& I1, Image<FVector<float, 3> >& I2, int iter_max, float theta, float alpha, float beta, float epsilon, float lambda)
{
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
	FVector<Image<FVector<float, 2>, 2 >, 3> gradI;
	gradI[0] = image_gradient(I1_R);
	gradI[1] = image_gradient(I1_G);
	gradI[2] = image_gradient(I1_B);

	// 2) Calcul de dtu
	// On calcule maintenant dtu pour les 3 composantes
	FVector<Image<float, 2>, 3> dtI;
	dtI[0] = image_dtu(I1, I2, 0);
	dtI[1] = image_dtu(I1, I2, 1);
	dtI[2] = image_dtu(I1, I2, 2);

	// 3) On initialise une carte de vitesse à 0
	FVector<Image<FVector<float, 2>, 2 >, 3> V;
	V[0] = init_map(w, h, 0.0f, 0.0f);
	V[1] = init_map(w, h, 0.0f, 0.0f);
	V[2] = init_map(w, h, 0.0f, 0.0f);

    // Vecteur de gradient de V
	FVector<Image<FVector<FVector<float, 2>, 2>, 2 >, 3> gradV;

	FVector<Image<FVector<float, 2>, 2 >, 3> W;
	W[0] = init_map(w, h);
	W[1] = init_map(w, h);
	W[2] = init_map(w, h);

	// Construisons le vecteur n a partir du gradient d'intensité
	FVector<Image<FMatrix<float, 2, 1>, 2>, 3> n;
	FVector<Image<FMatrix<float, 2, 1>, 2>, 3> n_ort;
	// Construisons le carré de la matrice D, nommé ici A
	FVector<Image<FMatrix<float, 2, 2>, 2>, 3> A;
	// Il reste a initialiser le vecteur p...
	FVector<Image<FVector<FVector<float, 2>, 2>, 2>, 3> p;
	FVector<Image<FVector<FVector<float, 2>, 2>, 2>, 3> B;

    
    // 4) On initialise les variables de boucle
	for (int comp = 0; comp < 3; comp++)
	{
        // allocation des images dans n, n_ort, A, B, p, et gradV
        n[comp] = Image<FMatrix<float, 2, 1>, 2>(w,h);
        n_ort[comp] = Image<FMatrix<float, 2, 1>, 2>(w,h);
        A[comp] = Image<FMatrix<float, 2, 2>, 2>(w,h);
        B[comp] = Image<FVector<FVector<float, 2>, 2>, 2>(w,h);
        p[comp] = Image<FVector<FVector<float, 2>, 2>, 2>(w,h);
        gradV[comp] = Image<FVector<FVector<float, 2>, 2>, 2 >(w,h);
        
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
                // allocation des matrices et vecteurs dans n, n_ort, A, B, p, et gradV
                n[comp](i, j) = FMatrix<float, 2, 1>();
                n_ort[comp](i, j) = FMatrix<float, 2, 1>();
                A[comp](i, j) = FMatrix<float, 2, 2>();
                B[comp](i, j) = FVector<FVector<float, 2>, 2>();

                p[comp](i, j) = FVector<FVector<float, 2>, 2>();
                p[comp](i, j)[0] = FVector<float, 2>();
                p[comp](i, j)[1] = FVector<float, 2>();

                gradV[comp](i, j) = FVector<FVector<float, 2>, 2>();
                gradV[comp](i, j)[0] = FVector<float, 2>();
                gradV[comp](i, j)[1] = FVector<float, 2>();
                
                float norm_gradI_ij = norm(gradI[comp](i,j));
				n[comp](i, j)(0,0) = (norm_gradI_ij > 0.0 )? gradI[comp](i, j)[0]/norm_gradI_ij : 0.0; // empêche la division par zero
				n[comp](i, j)(1,0) = (norm_gradI_ij > 0.0 )? gradI[comp](i, j)[1]/norm_gradI_ij : 0.0; // empêche la division par zero

				n_ort[comp](i, j)(0,0) = n[comp](i, j)[1];
				n_ort[comp](i, j)(1,0) = - n[comp](i, j)[0];

				// Represente le terme D^(1/2)
				A[comp](i, j) = exp(-alpha * pow(norm(gradI[comp](i, j)), beta))*(n[comp](i, j)*transpose(n[comp](i, j)) + n_ort[comp](i, j)*transpose(n_ort[comp](i, j)));
			}
		}
	}


	// 4 bis) Puis on met à jour l'estimation de la vitesse tant que le nombre d'iterations maximal n'est pas atteint
	
	// Utilisés dans les calculs :
	FVector<float,2> diff;
	float tau = 1.0 / (4.0 + epsilon);
	
	// iterations
	for (int iter = 0; iter < iter_max; iter++)
	{
		for (int comp = 0; comp < 3; comp++)
		{
			for (int i = 1; i < w - 1; i++)
			{
				for (int j = 1; j < h - 1; j++)
				{

					// Calcul de gradV
					gradV[comp](i, j) = gradient_2D(V[comp] ,i , j);

                    for (int dim = 0; dim < 2; dim++)
                    {
                    	// Actualisation de p
                    	diff = p[comp](i, j)[dim] * (1.0 - tau*epsilon)  +  tau * A[comp](i, j)*gradV[comp](i, j)[dim];
						float maxQ = max(1.0f, norm(diff) );
						cout << maxQ << " "<< norm(diff) << endl;
						p[comp](i, j)[dim] =  diff / maxQ;

						// Représente le produit D^(1/2)*p_{d}^{n+1}
						B[comp](i, j)[dim] = A[comp](i, j) * p[comp](i, j)[dim];

						// Actualisation de V
						V[comp](i, j)[dim] = W[comp](i, j)[dim] + theta * divergence(B[comp](i, j + 1)[dim], B[comp](i, j - 1)[dim], B[comp](i + 1, j)[dim], B[comp](i - 1, j)[dim]);
					}

					// Actualisation de W
					W[comp](i, j) = V[comp](i, j);

					float rho_ = rho(V[comp](i, j), gradI[comp](i, j), I1(i, j)[comp], I2(i, j)[comp]);
                	float norm_gradI_2 = gradI[comp](i, j) * gradI[comp](i, j);
                	float lambda_theta_normgradI_2 = rho_ - lambda * theta * norm_gradI_2 ;

					if (rho_ < -1.0 * lambda_theta_normgradI_2 )
						W[comp](i, j) += lambda * theta * gradI[comp](i, j);
					else if (rho_ > lambda_theta_normgradI_2 )
						W[comp](i, j) -= lambda * theta * gradI[comp](i, j);
					else{
						if( norm_gradI_2 != 0.0) // empêche la division par zero
							W[comp](i, j) -= rho_ * gradI[comp](i, j) / norm_gradI_2;
					}

					// Debug
					// cout << W[comp](i, j) << " ";
					// printf("%f\n",norm_gradI_2);

				}
			}
		}
        cout << "iteration : " << iter << " sur " << iter_max << endl;
	}

	// On retourne la moyenne de la vitesse sur chaque composante
	return 	(V[0] + V[1] + V[2])/3.0;;
}


// calcul un point fixe sur p et V
void point_fixe_V_p(FVector<float, 2>& V, FVector<FVector<float, 2>, 2>& p, int iter_max){
	// for (int i = 0; i < iter_max; i++){

	// }
}

// calcul la divergence d'un vecteur
float divergence(FVector<float, 2> b_droite, FVector<float, 2> b_gauche, FVector<float, 2> b_haut, FVector<float, 2> b_bas)
{
	return (b_bas[0] - b_haut[0])/2.0 + (b_droite[1] - b_gauche[1])/2.0;
}

// 
float rho(const FVector<float, 2>& u, const FVector<float, 2>& gradI, const float& I1, const float& I0)
{
	return u*gradI + I1 - I0;
}

// initialise une image de vecteurs aléatoirement
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

	for (int i=1; i<gradI.width()-1 ; i++){
		for (int j=1; j<gradI.height()-1 ; j++){  // pour chaque pixel, appel de Image::gradient
			gradI(i,j) = gradient_2D(I, i, j);
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

// Calcule la moyenne des vitesses autour du point i,j pour la composante c
FVector<float, 2> V_moy(Image<FVector<float, 2>, 2 >& V_c, int i, int j){
	return (2 * (V_c(i + 1, j) + V_c(i, j + 1) + V_c(i - 1, j) + V_c(i, j - 1)) + V_c(i + 1, j + 1) + V_c(i + 1, j - 1) + V_c(i - 1, j + 1) + V_c(i - 1, j - 1)) / 12.0;
}

// calcul du nouveau vecteur de vitesse en i,j pour la composante c
FVector<float, 2> nouveau_V_c( FVector<float, 2> V_moy_c_ij, FVector<float, 2> gradU_c_ij, float dtu_c_ij, float smoothness){
	return  V_moy_c_ij - gradU_c_ij * (gradU_c_ij * V_moy_c_ij + dtu_c_ij) / (smoothness*smoothness + gradU_c_ij * gradU_c_ij);
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

	// normalisation de la norme L2 des vecteurs entre 0 et 1
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
