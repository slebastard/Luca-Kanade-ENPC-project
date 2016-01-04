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

Image<FVector<float,2> , 2> find_flow_kanade_bw(Image<float>& I1, Image<float>& I2);
Image<FVector<FVector<float,2>,3> ,2 > find_flow_kanade(Image<FVector<float,3> >& I1, Image<FVector<float,3> >& I2);