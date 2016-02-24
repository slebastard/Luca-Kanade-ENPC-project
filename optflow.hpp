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
    double r;       // percent
    double g;       // percent
    double b;       // percent
} rgb;

typedef struct {
    double h;       // angle in degrees
    double s;       // percent
    double v;       // percent
} hsv;


Image<FVector<float,2> , 2> image_gradient(Image<float>& I);
Image<float, 2> image_dtu(Image<FVector<float,3> >& I1, Image<FVector<float,3> >& I2, int c);
Image<FVector<float,2> ,2 > optical_flow_calculation(Image<FVector<float,3> >& I1, Image<FVector<float,3> >& I2, int taille_fenetre, string method);
Image<FVector<float, 2> ,2 > flow_Lucas_Kanade(Image<FVector<float,3> >& I1, Image<FVector<float,3> >& I2, int taille_fenetre=7);
Image<FVector<float, 2>, 2 > flow_Horn_Schunk(Image<FVector<float, 3> >& I1, Image<FVector<float, 3> >& I2, float smoothness = 1, float stop = 0.01, int iter_max = 1000);
Image<FVector<float, 2>, 2 > flow_Horn_Schunk_HuberL1(Image<FVector<float, 3> >& I1, Image<FVector<float, 3> >& I2, int iter_max = 10000, float theta = 0.0002, float alpha = 1, float beta = 1, float epsilon = 0.0001);
Image<FVector<float, 2>, 2 > init_map(int width, int height, float v_min = 0, float v_max = 1);
FVector<float, 2> V_moy(Image<FVector<float, 2>, 2 >& V_c, int i, int j);
FVector<float, 2> nouveau_V_c( FVector<float, 2> V_moy_c_ij, FVector<float, 2> gradU_c_ij, float dtu_c_ij, float smoothness);
Image<Color, 2> make_flow_visible_hsv(Image<FVector<float,2> ,2 >& I);
Image<Color, 2> make_flow_visible_grey(Image<FVector<float,2> ,2 >& I);
rgb hsv2rgb(hsv in);
float angle_oriente(FVector<float,2>& v1, FVector<float,2>& v2);
float sgn(float val);

