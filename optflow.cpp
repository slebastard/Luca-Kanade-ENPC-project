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

#include "optflow.hpp"

using namespace std;
using namespace Imagine;


Image<FVector<FVector<float,2>,3> ,2 > find_flow_kanade(Image<FVector<float,3> >& I1, Image<FVector<float,3> >& I2){

	int h = I1.width(), w = I1.height();

	Image<FVector<FVector<float,2>,3> ,2 > V(w,h);
	Image<FVector<FVector<float,2>,3> ,2 > Gradu(w,h);
	Image<FVector<float,3> ,2 > dtu(w,h);


	// composante 1
	Image<float, 2> I1_R(w,h), I2_R(w,h);
	Image<float, 2> I1_G(w,h), I2_G(w,h);
	Image<float, 2> I1_B(w,h), I2_B(w,h);

	for (int i=0; i<I1.width() ; i++){
		for (int j=0; j<I1.height() ; j++){  // pour chaque pixel, appel de Image::gradient
			I1_R(i,j) = I1(i,j)[0];
			I1_G(i,j) = I1(i,j)[1];
			I1_B(i,j) = I1(i,j)[2];

			I2_R(i,j) = I2(i,j)[0];
			I2_G(i,j) = I2(i,j)[1];
			I2_B(i,j) = I2(i,j)[2];
		}
	}

	// calcul des composantes du gradient

	Image<FVector<float,2>, 2 > Gradu_R = find_flow_kanade_bw(I1_R, I2_R);
	Image<FVector<float,2>, 2 > Gradu_G = find_flow_kanade_bw(I1_G, I2_G);
	Image<FVector<float,2>, 2 > Gradu_B = find_flow_kanade_bw(I1_B, I2_B);

	// recomposer le gradient
	for (int i=0; i<Gradu.width() ; i++){
		for (int j=0; j<Gradu.height() ; j++){  // pour chaque pixel, appel de Image::gradient
			Gradu(i,j)[0] = Gradu_R(i,j);
			Gradu(i,j)[1] = Gradu_G(i,j);
			Gradu(i,j)[2] = Gradu_B(i,j);
		}
	}

	// On a maintenant le gradient


}

Image<FVector<float,2> , 2> find_flow_kanade_bw(Image<float>& I1, Image<float>& I2){

	int h = I1.width(), w = I1.height();

	Image<FVector<float,2> ,2> V(w,h);
	Image<FVector<float,2> ,2> Gradu(w,h);
	Image<float, 2> dtu(w,h);

	// calcul de Gradu

	for (int i=0; i<w ; i++){
		for (int j=0; j<h ; j++){  // pour chaque pixel, appel de Image::gradient
			Gradu(i,j) = gradient( I1, Coords<2>(i,j) );
		}
	}
	
}
