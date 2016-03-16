#ifndef OPTFLOW_H
#define OPTFLOW_H

/*
Projet MOPSI : Flux optique
Loic Cressot  &  Simon Lebastard
Encadrant : Pascal Monasse
janvier 2016
*/

#include <iostream>
#include <string>
#include <vector>

#include <Imagine/Images.h>
#include "Imagine/Common.h"

using namespace std;
using namespace Imagine;


typedef struct {
    double r;       // pourcentage
    double g;       // pourcentage
    double b;       // pourcentage
} rgb;

typedef struct {
    double h;       // angle en degres
    double s;       // pourcentage
    double v;       // pourcentage
} hsv;


Image<double, 2> image_dtu(Image<double >& I1, Image<double >& I2);
Image<double, 2> image_dtu(Image<FVector<double,3> >& I1, Image<FVector<double,3> >& I2, int c);

Image<FVector<double, 2> ,2 > flow_Lucas_Kanade(Image<FVector<double,3> >& I1, Image<FVector<double,3> >& I2, int taille_fenetre=7);
Image<FVector<double, 2>, 2 > flow_Horn_Schunk(Image<FVector<double, 3> >& I1, Image<FVector<double, 3> >& I2, double smoothness = 1.0, double stop = 0.01, int iter_max = 1000, bool blackNwhite=false);
Image<FVector<double, 2>, 2 > flow_Horn_Schunk_Grey(Image<double>& I1, Image<double>& I2, double smoothness, double stop, int iter_max);
Image<FVector<double, 2>, 2 > flow_Horn_Schunk_HuberL1(Image<FVector<double, 3> >& I1, Image<FVector<double, 3> >& I2, int iter_max = 50, double lambda = 1.0, double theta = 0.0002, double alpha = 1.0, double beta = 1.0, double epsilon = 0.1);

Image<FVector<double, 2>, 2 > init_map(int width, int height, double v_min = 0, double v_max = 1);
double divergence(FVector<double, 2> b_droite, FVector<double, 2> b_gauche, FVector<double, 2> b_haut, FVector<double, 2> b_bas);
double rho(const FVector<double, 2>& u, const FVector<double, 2>& gradI, const double& I1, const double& I0);
FVector<double, 2> V_moy(Image<FVector<double, 2>, 2 >& V_c, int i, int j);
FVector<double, 2> nouveau_V_c( FVector<double, 2> V_moy_c_ij, FVector<double, 2> gradU_c_ij, double dtu_c_ij, double smoothness);

void convert_to_BW(Image<FVector<double, 3> >& I);
Image<Color, 2> make_flow_visible_hsv(Image<FVector<double,2> ,2 >& I);
Image<Color, 2> make_flow_visible_grey(Image<FVector<double,2> ,2 >& I);

rgb hsv2rgb(hsv in);
double angle_oriente(FVector<double,2>& v1, FVector<double,2>& v2);
double sgn(double val);
long double normGradV(Image<FVector<FVector<double, 2>, 2> ,2 >& gV);


// calcule le gradient en i,j d'une image de T (gradient 2D approximé à l'ordre 1 : cf filtre de Sobel)
// prendre en compte les bords pour amélioration
template <typename T>
FVector<T, 2> gradient_2D(const Image<T, 2> &V, int i, int j) {
	FVector<T, 2> gradV;
	// composante x
	gradV[0] = (   ( V(i-1,j-1) + 2*V(i-1,j) + V(i-1,j+1) )*(-1.0)
			     + ( V(i+1,j-1) + 2*V(i+1,j) + V(i+1,j+1) )
			   ) ;
	
	// composante y
	gradV[1] = (   ( V(i-1,j-1) + 2*V(i,j-1) + V(i+1,j-1) )*(-1.0)
			     + ( V(i-1,j+1) + 2*V(i,j+1) + V(i+1,j+1) )
			   ) ;

	return gradV;
}

// Calcule le gradient 2D pour toute une image
template <typename T>
Image<FVector<T, 2>, 2> image_gradient_2D(const Image<T, 2> &V) {

	Image<FVector<T, 2>, 2> gradV(V.width(),V.height());
	
	for(int i=1; i<V.width()-1; i++){
		for(int j=1; j<V.height()-1; j++){
			gradV(i,j) = gradient_2D(V, i, j);
		}
	}

	return gradV;
}


// Calcule un produit scalaire de deux images template
template <typename T>
Image<T, 2> image_dot(const Image<FVector<T, 2>, 2> &I1, const Image<FVector<T, 2>, 2> &I2){
	Image<T, 2> dotp(I1.width(),I1.height());
	for(int i=1; i<I1.width()-1; i++){
		for(int j=1; j<I1.height()-1; j++){
			dotp(i,j) = I1(i,j)*I2(i,j);
		}
	}
    return dotp;
}

#endif