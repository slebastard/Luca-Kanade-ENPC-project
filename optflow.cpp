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
#include <Imagine/Common.h>
#include <Imagine/LinAlg.h>

#include "optflow.hpp"

using namespace std;
using namespace Imagine;

// Calcul du flow optique par la methode Lucas et Kanade
Image<FVector<double,2> ,2 > flow_Lucas_Kanade(Image<FVector<double,3> >& I1, Image<FVector<double,3> >& I2, int taille_fenetre){

	int w = I1.width(), h = I1.height();
	int w2 = I2.width(), h2 = I2.height();

	if(w!=w2 || h!=h2){
		// Erreur si images de taille différentes
		throw string("Erreur : Images de tailles différentes");
	}

	// 1) Calcul du gradient
	// On decompose l'image 1 en trois images selon ses 3 composantes R G B 
	Image<double, 2> I1_R(w,h);
	Image<double, 2> I1_G(w,h);
	Image<double, 2> I1_B(w,h);

	for (int i=0; i < w ; i++){
		for (int j=0; j < h ; j++){ 
			I1_R(i,j) = I1(i,j)[0];
			I1_G(i,j) = I1(i,j)[1];
			I1_B(i,j) = I1(i,j)[2];
		}
	}

	// Calcul des composantes des images gradient des 3 composantes R G B
	Image<FVector<double,2>, 2 > gradU_R = image_gradient_2D(I1_R);
	Image<FVector<double,2>, 2 > gradU_G = image_gradient_2D(I1_G);
	Image<FVector<double,2>, 2 > gradU_B = image_gradient_2D(I1_B);

	// 2) Calcul de dtu
	// On calcule maintenant dtu pour les 3 composantes
	Image<double, 2> dtu_R = image_dtu(I1,I2,0);
	Image<double, 2> dtu_G = image_dtu(I1,I2,1);
	Image<double, 2> dtu_B = image_dtu(I1,I2,2);

	//3) Caclcul du flow optique par régression linéaire
	//On crée l'image qui contiendra le flow optique que l'on va retourner
	Image<FVector<double,2> ,2 > V(w,h); // sera retourné une fois rempli

	// Calcul flow optique :
	int n = (int)taille_fenetre/2;
    const int n_equations = taille_fenetre * taille_fenetre * 3;

	Matrix<double> A(n_equations,2);
	Vector<double> b(n_equations);

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
            Vector<double> v = linSolve(A,b);
            V(i,j)[0] = v[0];
            V(i,j)[1] = v[1];
		}
	}

	// On retourne le flow optique
	return V;
}



// Calcul du flow optique par la methode Horn et Schunk itérative pour des images en couleur
Image<FVector<double, 2>, 2 > flow_Horn_Schunk(Image<FVector<double, 3> >& I1, Image<FVector<double, 3> >& I2, double smoothness, double stop, int iter_max, bool blackNwhite){

	// COnvertit les image en nuances de gris (si demandé)
	if(blackNwhite){
		convert_to_BW(I1);
		convert_to_BW(I2);
	}

	int w = I1.width(), h = I1.height();
	int w2 = I2.width(), h2 = I2.height();

	if (w != w2 || h != h2){
		// Erreur si images de taille différentes
		throw string("Erreur : Images de tailles différentes");
	}

	// 1) On decompose les images en trois images selon les 3 composantes R G B 
	Image<double, 2> I1_R(w, h);
	Image<double, 2> I1_G(w, h);
	Image<double, 2> I1_B(w, h);

	Image<double, 2> I2_R(w, h);
	Image<double, 2> I2_G(w, h);
	Image<double, 2> I2_B(w, h);

	for (int i = 0; i < w; i++){
		for (int j = 0; j < h; j++){
			I1_R(i, j) = I1(i, j)[0];
			I1_G(i, j) = I1(i, j)[1];
			I1_B(i, j) = I1(i, j)[2];
			I2_R(i, j) = I2(i, j)[0];
			I2_G(i, j) = I2(i, j)[1];
			I2_B(i, j) = I2(i, j)[2];
		}
	}

	// Et on applique HS sur chaque composantes
	Image<FVector<double, 2>, 2 > V_R = flow_Horn_Schunk_Grey(I1_R, I2_R, smoothness, stop, iter_max);
	Image<FVector<double, 2>, 2 > V_G = flow_Horn_Schunk_Grey(I1_G, I2_G, smoothness, stop, iter_max);
	Image<FVector<double, 2>, 2 > V_B = flow_Horn_Schunk_Grey(I1_B, I2_B, smoothness, stop, iter_max);

	//3.3) On renvoie la moyenne des 3 résultats
	return  ( V_R + V_G + V_B ) / 3.0 ;
}


// Calcul du flow optique par la methode Horn et Schunk itérative sur des images en nuances de gris
Image<FVector<double, 2>, 2 > flow_Horn_Schunk_Grey(Image<double>& I1, Image<double>& I2, double smoothness, double stop, int iter_max){

	
	int w = I1.width(), h = I1.height();
	int w2 = I2.width(), h2 = I2.height();

	if (w != w2 || h != h2){
		// Erreur si images de taille différentes
		throw string("Erreur : Images de tailles différentes");
	}

	// 1) Calcul du gradient

	// Calcul des composantes des images gradient
	Image<FVector<double, 2>, 2 > gradU = image_gradient_2D(I1);

	// 2) Calcul de dtu
	// On calcule maintenant dtu
	Image<double, 2> dtu = image_dtu(I1, I2);

	/*   Calcul du flow optique par itération sur le modèle de Horn & Schunk    */

	//3.1) Créer une carte de vitesses initiales, idéalement en prenant comme entrée le résultat d'une autre méthode
	//     Pour l'instant l'entrée en vecteurs vitesses sera aléatoire
	Image<FVector<double, 2>, 2 > V = init_map(w, h);
	Image<FVector<double, 2>, 2 > V_ant(w, h);
	Image<FVector<double, 2>, 2 > V_m(w, h);


	int iter = 0;
	double energie, energie_ant;
	//3.2) 
	while (iter < iter_max)
	{
		energie_ant = energie;
		for (int i = 1; i < w-1; i++){
			for (int j = 1; j < h-1; j++){
				// calcul de la moyenne locale R
				V_m(i, j) = V_moy(V,i,j);
				// calcul du nouveau V
				V(i,j) = nouveau_V_c( V_m(i, j), gradU(i, j), dtu(i,j), smoothness);
			}
		}
		// on affiche l'énergie :
        Image<FVector<FVector<double, 2>, 2>, 2 > gradV = image_gradient_2D(V);
        energie = norm( image_dot(gradU,V) + dtu) + normGradV(gradV);
		//cout << iter << " " << std::to_string(energie) << endl;
		iter++;
		// condition d'arrêt supplémentaire : convergence de l'énergie
		if( abs(energie_ant - energie) < stop) break;
	}
    cout << iter-1 << " itérations" << endl;
	
	// On retourne le flow optique
	return V;
}


// Calcul du flow optique par la methode Horn et Schunk Huber L1 (non encore fonctionnelle)
Image<FVector<double, 2>, 2 > flow_Horn_Schunk_HuberL1(Image<FVector<double, 3> >& I1, Image<FVector<double, 3> >& I2, int iter_max, double theta, double alpha, double beta, double epsilon, double lambda)
{    
	int w = I1.width(), h = I1.height();
	int w2 = I2.width(), h2 = I2.height();

	if (w != w2 || h != h2){
		// Erreur si images de taille différentes
		throw string("Erreur : Images de tailles différentes");
	}

	// 1) Calcul du gradient
	// On decompose l'image 1 en trois images selon ses 3 composantes R G B 
	Image<double, 2> I1_R(w, h);
	Image<double, 2> I1_G(w, h);
	Image<double, 2> I1_B(w, h);

	for (int i = 0; i < w; i++){
		for (int j = 0; j < h; j++){
			I1_R(i, j) = I1(i, j)[0];
			I1_G(i, j) = I1(i, j)[1];
			I1_B(i, j) = I1(i, j)[2];
		}
	}

	// Calcul des composantes des images gradient des 3 composantes R G B
	FVector<Image<FVector<double, 2>, 2 >, 3> gradI;
	gradI[0] = image_gradient_2D(I1_R);
	gradI[1] = image_gradient_2D(I1_G);
	gradI[2] = image_gradient_2D(I1_B);

	// 2) Calcul de dtu
	// On calcule maintenant dtu pour les 3 composantes
	FVector<Image<double, 2>, 3> dtI;
	dtI[0] = image_dtu(I1, I2, 0);
	dtI[1] = image_dtu(I1, I2, 1);
	dtI[2] = image_dtu(I1, I2, 2);

	// 3) On initialise une carte de vitesse à 0
	FVector<Image<FVector<double, 2>, 2 >, 3> V;
	V[0] = init_map(w, h, 0.0f, 0.0f);
	V[1] = init_map(w, h, 0.0f, 0.0f);
	V[2] = init_map(w, h, 0.0f, 0.0f);

    // Vecteur de gradient de V
	FVector<Image<FVector<FVector<double, 2>, 2>, 2 >, 3> gradV;

	FVector<Image<FVector<double, 2>, 2 >, 3> W;
	W[0] = init_map(w, h);
	W[1] = init_map(w, h);
	W[2] = init_map(w, h);

	// Construisons le vecteur n a partir du gradient d'intensité
	FVector<Image<FMatrix<double, 2, 1>, 2>, 3> n;
	FVector<Image<FMatrix<double, 2, 1>, 2>, 3> n_ort;
	// Construisons le carré de la matrice D, nommé ici A
	FVector<Image<FMatrix<double, 2, 2>, 2>, 3> A;
	// Il reste a initialiser le vecteur p...
	FVector<Image<FVector<FVector<double, 2>, 2>, 2>, 3> p;
	FVector<Image<FVector<FVector<double, 2>, 2>, 2>, 3> B;

    
    // 4) On initialise les variables de boucle
	for (int comp = 0; comp < 3; comp++)
	{
        // allocation des images dans n, n_ort, A, B, p, et gradV
        n[comp] = Image<FMatrix<double, 2, 1>, 2>(w,h);
        n_ort[comp] = Image<FMatrix<double, 2, 1>, 2>(w,h);
        A[comp] = Image<FMatrix<double, 2, 2>, 2>(w,h);
        B[comp] = Image<FVector<FVector<double, 2>, 2>, 2>(w,h);
        p[comp] = Image<FVector<FVector<double, 2>, 2>, 2>(w,h);
        gradV[comp] = Image<FVector<FVector<double, 2>, 2>, 2 >(w,h);
        
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
                // allocation des matrices et vecteurs dans n, n_ort, A, B, p, et gradV
                n[comp](i, j) = FMatrix<double, 2, 1>();
                n_ort[comp](i, j) = FMatrix<double, 2, 1>();
                A[comp](i, j) = FMatrix<double, 2, 2>();
                B[comp](i, j) = FVector<FVector<double, 2>, 2>();

                p[comp](i, j) = FVector<FVector<double, 2>, 2>();
                p[comp](i, j)[0] = FVector<double, 2>();
                p[comp](i, j)[1] = FVector<double, 2>();

                gradV[comp](i, j) = FVector<FVector<double, 2>, 2>();
                gradV[comp](i, j)[0] = FVector<double, 2>();
                gradV[comp](i, j)[1] = FVector<double, 2>();
                
                double norm_gradI_ij = norm(gradI[comp](i,j));
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
	FVector<double,2> diff;
	double tau = 1.0 / (4.0 + epsilon);
	
	// iterations
	for (int iter = 0; iter < iter_max; iter++)
    { cout << endl << endl;
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
                        double maxQ = std::max(1.0, norm(diff) );
						//cout << maxQ << " "<< norm(diff) << endl;
						p[comp](i, j)[dim] =  diff / maxQ;

						// Représente le produit D^(1/2)*p_{d}^{n+1}
						B[comp](i, j)[dim] = A[comp](i, j) * p[comp](i, j)[dim];
                        
						// Actualisation de V
						double a = W[comp](i, j)[dim] + theta * divergence(B[comp](i, j + 1)[dim], B[comp](i, j - 1)[dim], B[comp](i + 1, j)[dim], B[comp](i - 1, j)[dim]);
                        V[comp](i, j)[dim] = a;
                        if(abs(a)>1000){ cout << "a iter i j"<< a << " " << iter << " " <<  i <<  " " << j << " " <<  endl;
                            exit(1);
                        }
					}

					// Actualisation de W

					double rho_ = rho(V[comp](i, j), gradI[comp](i, j), I1(i, j)[comp], I2(i, j)[comp]);
                	double norm_gradI_2 = gradI[comp](i, j) * gradI[comp](i, j);
                	double lambda_theta_normgradI_2 = rho_ - lambda * theta * norm_gradI_2 ;


					W[comp](i, j) = V[comp](i, j);
					if (rho_ < -1.0 * lambda_theta_normgradI_2 )
						W[comp](i, j) += lambda * theta * gradI[comp](i, j);
					else if (rho_ > lambda_theta_normgradI_2 )
						W[comp](i, j) -= lambda * theta * gradI[comp](i, j);
					else{
						if( norm_gradI_2 != 0.0) // empêche la division par zero
							W[comp](i, j) -= rho_ * gradI[comp](i, j) / norm_gradI_2;
					}
				}
			}
		}
        cout << "iteration : " << iter << endl;
    }
    
	// On retourne la moyenne de la vitesse sur chaque composante
	return 	(V[0] + V[1] + V[2])/3.0;;
}





// Calcul la divergence d'un vecteur
double divergence(FVector<double, 2> b_droite, FVector<double, 2> b_gauche, FVector<double, 2> b_haut, FVector<double, 2> b_bas)
{
	return (b_bas[0] - b_haut[0])/2.0 + (b_droite[1] - b_gauche[1])/2.0;
}

// Calcul de rho
double rho(const FVector<double, 2>& u, const FVector<double, 2>& gradI, const double& I1, const double& I0)
{
	return u*gradI + I1 - I0;
}

// Initialise une image de vecteurs aléatoirement
Image<FVector<double, 2>, 2 > init_map(int width, int height, double v_min, double v_max)
{
	Image<FVector<double, 2>, 2 > V(width, height);
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


// Calcul de dtu la dérivée temporelle entre deux images en nuances de gris
Image<double, 2> image_dtu(Image<double >& I1, Image<double >& I2){
	
	Image<double, 2> dtu(I1.width(),I1.height());
	
	for (int i=0; i<dtu.width() ; i++){
		for (int j=0; j<dtu.height() ; j++){  // pour chaque pixel, appel de Image::gradient
			dtu(i,j) = I2(i,j) - I1(i,j);
		}
	}
	return dtu;
}

// Calcul de dtu la dérivée temporelle pour une composante donnée entre deux images couleurs
Image<double, 2> image_dtu(Image<FVector<double,3> >& I1, Image<FVector<double,3> >& I2, int c){
	
	Image<double, 2> dtu(I1.width(),I1.height());
	
	for (int i=0; i<dtu.width() ; i++){
		for (int j=0; j<dtu.height() ; j++){  // pour chaque pixel, appel de Image::gradient
			dtu(i,j) = I2(i,j)[c] - I1(i,j)[c];
		}
	}
	return dtu;
}

// Calcule la moyenne des vitesses autour du point i,j pour la composante c
FVector<double, 2> V_moy(Image<FVector<double, 2>, 2 >& V_c, int i, int j){
	return (2 * (V_c(i + 1, j) + V_c(i, j + 1) + V_c(i - 1, j) + V_c(i, j - 1)) + V_c(i + 1, j + 1) + V_c(i + 1, j - 1) + V_c(i - 1, j + 1) + V_c(i - 1, j - 1)) / 12.0;
}

// calcul du nouveau vecteur de vitesse en i,j pour la composante c
FVector<double, 2> nouveau_V_c( FVector<double, 2> V_moy_c_ij, FVector<double, 2> gradU_c_ij, double dtu_c_ij, double smoothness){
	return  V_moy_c_ij - gradU_c_ij * (gradU_c_ij * V_moy_c_ij + dtu_c_ij) / (smoothness*smoothness + gradU_c_ij * gradU_c_ij);
}

// Converti une image couleur en nuances de gris (image couleur avec composantes identiques et grises)
void convert_to_BW(Image<FVector<double, 3> >& I){
	for(int i=0; i<I.width(); i++){
		for(int j=0; j<I.height(); j++){
			I(i,j)[0] = I(i,j)[1] = I(i,j)[2] = 0.21f * I(i,j)[0] + 0.72f * I(i,j)[1] + 0.07f * I(i,j)[2];
		}
	}
}


// Transforme le flow optique en image visualisable
Image<Color, 2> make_flow_visible_hsv(Image<FVector<double,2> ,2 >& I){

	int w = I.width(), h = I.height();
	Image<Color, 2> img(w,h);

	// recherche de la norme maximum
	double maximum = -1.0*FLT_MAX;
	for (int i=0; i<w ; i++){
		for (int j=0; j<h ; j++){
			maximum = max(maximum, I(i,j)*I(i,j));
		}
	}

	// normalisation de la norme L2 des vecteurs entre 0 et 1
	I/=sqrt(maximum);

	// reference vector :
	FVector<double,2> ref_vect(0.0,1.0);
	// Maintenant on peut donner une teinte et valeur à chaque pixel de img avec les valeurs de I
	for (int i=0; i<w ; i++){
		for (int j=0; j<h ; j++){  // pour chaque pixel, appel de Image::gradient
			hsv in;
			double norm = sqrt(I(i,j)*I(i,j)); // Okay
			double angle = angle_oriente(I(i,j),ref_vect)*180.0/M_PI + 180.0; // okay
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
Image<Color, 2> make_flow_visible_grey(Image<FVector<double,2> ,2 >& I){

	int w = I.width(), h = I.height();
	Image<Color, 2> img(w,h);

	double maximum = -1.0*FLT_MAX;
	double minimum = FLT_MAX;
	for (int i=0; i<w ; i++){
		for (int j=0; j<h ; j++){  // pour chaque pixel, appel de Image::gradient
			maximum = max(maximum, I(i,j)[0]);
			maximum = max(maximum, I(i,j)[1]);
			minimum = min(minimum, I(i,j)[0]);
			minimum = min(minimum, I(i,j)[1]);
		}
	}

	maximum-= minimum;
	I-= FVector<double,2>(minimum,minimum);
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


// conversion hsv vers rgb
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
double angle_oriente(FVector<double,2>& v1, FVector<double,2>& v2){
	double L = doubleNorm(v1) * doubleNorm(v2);
	if(L < FLT_EPSILON*5)
		return 0.0;
	else{
		double angle = acos(v1*v2/L);
		double signe = sgn(v1^v2);
		return angle * signe;
	}
}

// Fonction signe
double sgn(double val){
	if(val < 0.0)
    	return -1.0;
    else
    	return 1.0;
}

// Calcule la norme d'un gradient de flow
long double normGradV(Image<FVector<FVector<double, 2>, 2> ,2 >& gV){
	 long double normG = 0.0;
	 for (int i=0; i<gV.width() ; i++){
	 	for (int j=0; j<gV.height() ; j++){  // pour chaque pixel, appel de Image::gradient
	 		normG += sqrt( gV(i,j)[0]*gV(i,j)[0] + gV(i,j)[1]*gV(i,j)[1] );
	 	}
	 }
    return normG;
}

