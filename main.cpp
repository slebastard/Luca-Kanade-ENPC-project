/*
Projet MOPSI : Flux optique
Loic Cressot  &  Simon Lebastard
Encadrant : Pascal Monasse
Janvier 2016
*/

#define S_IRWXU (S_IRUSR | S_IWUSR | S_IXUSR)
#define S_IRWXG (S_IRWXU >> 3)

#define NOMINMAX

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
  #include <windows.h>
  #include "getopt.h"
  #include "unistd.h"
  #include "dirent.h"
#else
  #include <getopt.h>
  #include <unistd.h>
  #include <dirent.h>
#endif


#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <sys/stat.h>
#include <sys/types.h>

#include <Imagine/Images.h>
#include "optflow.hpp"
#include "readflow.h"


using namespace std;
using namespace Imagine;


void showhelpinfo(string s = "optflow");


/*
  QUESTIONS :
  - paramètre int à un template
  - grad spatial pour chaque composantes ou 3d ?
  - moyenne sur les composantes ?
*/


int main (int argc,char *argv[])
{
    string directory("");
    string output_directory("");
	string ground_truth_path("");
    vector<string> method_args;
    bool verbose = false;
    bool save_outputs = false;
    bool print_outputs = false;
    bool gif_style = false;
	bool comp_error_map = false;
    int MAX_RES = 100000;


// PROCESS ARGUMENTS
// ===================================================
    char tmp;
/*if the program is ran without options ,it will show the usgage and exit*/
    if(argc == 1)
    {
        showhelpinfo(string(argv[0]));
        return EXIT_FAILURE;
    }
// get options
  while ((tmp = getopt(argc, argv, "hvsgpd:r:o:m:")) != -1)
    {
      switch (tmp)
      {
      /*option h show the help infomation*/
      case 'h':
        showhelpinfo(string(argv[0]));
        break;
      /*option v for verbose*/
      case 'v':
        verbose = true;
        break;
      /*option s for save*/
      case 's':
        save_outputs = true;
        break;
      /*option s for save*/
      case 'r':
        MAX_RES = atoi(optarg);
        break;
      /*option o for output dir*/
      case 'o':
        output_directory = string(optarg);
        break;
      /*option g for gif style format*/
      case 'g':
        gif_style = true;
        break;
      /*option p for print_outputs*/
      case 'p':
        print_outputs = true;
        break;
      /*option d asks for directory*/
      case 'd':
        directory = string(optarg);
        break;
	  /*option e for error estimation output*/
	  case 'e':
		ground_truth_path = string(optarg);
		comp_error_map = true;
		break;
      /*option m asks for method args*/
      case 'm':
        {
        istringstream iss(optarg);
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(method_args));
        break;
        }
      // default shows help
      default:
        showhelpinfo(string(argv[0]));
        break;
      }
    }

// verify options

//directory
if(directory==""){
    cout << "Directory option required !" << endl;
    return EXIT_FAILURE;
}
if(output_directory==""){
    output_directory = directory + "/outputs";
}
// save outputs
if(!save_outputs) print_outputs = true;
// method
if(method_args.size()==0){
    cout << "Method option required !" << endl;
    return EXIT_FAILURE;
}

if (method_args[0] != "LK" && method_args[0] != "HS" && method_args[0] != "HSL1"){
    cout << "Unknown method : " << method_args[0] << endl;
    showhelpinfo();
    return EXIT_FAILURE;
}
if(method_args[0]=="LK"){
  if(method_args.size()!=2){
    cout << "Incorrect args number for LK method ! " << endl;
    showhelpinfo();
    return EXIT_FAILURE;
  }
}
else if(method_args[0]=="HS"){
  if(method_args.size()!=4){
    cout << "Incorrect args number for HS method ! " << endl;
    showhelpinfo();
    return EXIT_FAILURE;
  }
}
else if (method_args[0] == "HSL1"){
	if (method_args.size() != 2){
		cout << "Incorrect args number for HSL1 method ! " << endl;
    showhelpinfo();
		return EXIT_FAILURE;
	}
}


// CHARGEMENT DES IMAGES
// ===================================================
cout << endl;
cout << "==================================================="<<endl;
cout << "CHARGEMENT DES IMAGES" << endl;
cout << "==================================================="<<endl;
// vector of images
vector<Image<FVector<float,3>, 2 > > images;
vector<Image<Color, 2 > > outputs;

//open directory and load images
DIR *dir;
struct dirent *ent;
if ((dir = opendir (directory.c_str())) != NULL) {
/* print all the files and directories within directory */
    while ((ent = readdir (dir)) != NULL) {
        //printf ("%s\n", ent->d_name);
        string name(ent->d_name);
        // look for supported formats
        if(name.find(".png")!=string::npos || name.find(".jpg")!=string::npos || name.find(".tiff")!=string::npos ){
            Image<Color> I;
            if(!load(I, directory + "/" + name)){
                cout << string("Couldn't load image : ") + name << endl;
                return EXIT_FAILURE;
            }
            Image<FVector<float,3> > Iv(I);
            images.push_back(Iv);
            if(verbose) cout << "Image correctly loaded : " << name << endl;
        }
        else{
            if(verbose) cout << "Ignoring uncompatible format file : " << name << endl;
        }
    }
    closedir (dir);
} else {
/* could not open directory */
    cout << string("Could not open directory ") + directory << endl;
    return EXIT_FAILURE;
}

// On vérifie le nombre d'images
if(images.size()<2){
    cout << "Not enough images loaded, need at least 2" << endl;
    return EXIT_FAILURE;
}



// PRETRAITEMENT DES IMAGES
// ===================================================
cout << endl;
cout << "==================================================="<<endl;
cout << "PRETRAITEMENT DES IMAGES" << endl;
cout << "==================================================="<<endl;

int w = images[0].width(), h = images[0].height();

// On vérifie que les images ont la même taille et ne sont pas trop grandes auquel cas on les réduit
for(int i=0; i<images.size(); i++){
  if(w != images[i].width() || h != images[i].height()){
    // Erreur si images de taille différentes
    throw string("Erreur : Images de tailles différentes");
  }
  if(images[i].width() * images[i].height() > MAX_RES){
    // Rescale si on depasse la resolution maximale
    double fact = 1.0*MAX_RES / (images[i].width() * images[i].height());
    if(verbose) cout << "Reducing image n°" << i << endl;
    images[i] = reduce(images[i], 1/fact);
  }

}







// TRAITEMENT DES IMAGES
// ===================================================
cout << endl;
cout << "==================================================="<<endl;
cout << "TRAITEMENT DES IMAGES ..." << endl;
cout << "==================================================="<<endl;


bool first = true;
for(int i=0; i<images.size()-1; i++){
  // Calcul du flow optique
  if(verbose) cout << "Processing image n°" << i << endl;

    Image<FVector<float,2> ,2 > optical_flow;
    
	if (method_args[0] == "LK")
	{
		int taille_fenetre = atoi(method_args[1].c_str());
		optical_flow = flow_Lucas_Kanade(images[i], images[i + 1], 7);
	}
	else if (method_args[0] == "HS")
	{
		float smoothness = atof(method_args[1].c_str());
		float stop = atof(method_args[2].c_str());
		int max_iter = atoi(method_args[3].c_str());
		optical_flow = flow_Horn_Schunk(images[i], images[i + 1], smoothness, stop, max_iter);
	}
	else if (method_args[0] == "HSL1")
	{
		int max_iter = atoi(method_args[1].c_str());
		optical_flow = flow_Horn_Schunk_HuberL1(images[i], images[i + 1], max_iter);
	}
    

  // Visualisation
  Image<Color, 2 > optical_flow_image = make_flow_visible_hsv(optical_flow);
  outputs.push_back(optical_flow_image);
  
  if(first && print_outputs){
    openWindow(optical_flow_image.width(), optical_flow_image.height()); 
    first=false;
  }
  if(print_outputs){
    display(optical_flow_image);
    anyClick();
  }
}







// SAUVEGARDE DES RESULTATS
// ===================================================
if(save_outputs){
  cout << endl;
  cout << "==================================================="<<endl;
  cout << "SAUVEGARDE DES RESULTATS" << endl;
  cout << "==================================================="<<endl;
}

if (save_outputs){

  // ouverture ou creation du dossier outputs

  if ((dir = opendir(output_directory.c_str())) == NULL)
  {
    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
      if (int e = CreateDirectory(output_directory.c_str(), NULL) != 0){
        throw string("Cannot make directory " + output_directory);
      }
    #else
      if ( mkdir( output_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
        throw string("Cannot make directory " + output_directory);
      }
    #endif
  }

  if(verbose) cout << "Saving images to " + output_directory << endl;
  int n_dec = int(log10(outputs.size()) + 1); // number of decimal to use for number in name
  char* num = new char[n_dec]; // tableau contenant le numéro  de l'image en chaine de caractères
  for(int i=0; i<outputs.size(); i++){
    cout << "saving " << "/output_" << to_string(i) << ".jpg" << endl;
      
    sprintf(num, (string("%0"+to_string(n_dec)+"d")).c_str(), i);
    save(outputs[i], output_directory + "/output_" + num +".jpg");
  }
  if(gif_style){
    int j=outputs.size();
    for(int i=outputs.size()-1; i>=0; i--){ 
      cout << "saving " << "output_" << to_string(j) << ".png" << endl;
      save(outputs[i], output_directory + "/output_" + to_string(j)+".png");
      j++;
    }
  }


}

// GENERATION DE LA CARTE D'ERREUR
// ===================================================
if (comp_error_map){
	// A COMPLETER
	// GENERATION DE LA CARTE D'ERREUR 
	// + TRAITEMENT ET AFFICHAGE
}


return 0;
}


/*fonction that show the help information*/
void showhelpinfo(string s)
{
    cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
    cout<<"option:  "<<"-h  show help information"<<endl;
    cout<<"         "<<"-d directory"<<endl;
    cout<<"         "<<"-o output directory"<<endl;
    cout<<"         "<<"-s save"<<endl;
    cout<<"         "<<"-p print results"<<endl;
    cout<<"         "<<"-m method args \"LK taille_fenetre\"  or\n           \"HS smoothness stop max_iter\"  or\n        \"HSL1 max_iter\" "<<endl;
    cout<<"         "<<"-v verbose"<<endl;
    cout<<"         "<<"-r verbose"<<endl;
    cout<<"         "<<"-r max resolution of ouptputs"<<endl;
    cout<<"         "<<"-g gif_style"<<endl;
    cout<<"example: "<<s<<" -d directory -s1"<<endl;
}
