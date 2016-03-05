#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

#include <Imagine/Images.h>
#include "Imagine/Common.h"
using namespace Imagine;

int hex2dec(string hex_byte);
int get_bytes(istream& stream, int nb_byte = 4);
bool right_stream_pos(istream& stream, int row, int col, int width, int height);

Image<FVector<float, 2>, 2 > flow_from_flo(string& path);

Image<float, 2> dim_wise_error(Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation, int dim);
Image<float, 2> norm_error(Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation, bool quadratic = false);
Image<float, 2> error_map(string feature, Image<FVector<float, 2>, 2> ground_truth, Image<FVector<float, 2>, 2> flow_estimation);
