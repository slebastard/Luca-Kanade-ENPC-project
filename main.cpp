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
#include <Imagine/Common.h>

#include "optflow.hpp"
#include "readflow.hpp"
#include "flowcode/flowIO.h"


using namespace std;
using namespace Imagine;


void showhelpinfo(string s = "optflow");
void write_lines_to_csv(const string& file_path, vector<string>& lines);



int main (int argc,char *argv[])
{
    // Variables d'options
    string directory("");         // stores input dir
    string output_directory("");  // stores output dir
    string ground_truth_path(""); // stores path to ground_truth file
    vector<string> method_args;   // stores args provided for options -m (method) between ""
    bool verbose = false;         // verbose
    bool save_outputs = false;    // save image outputs to output dir
    bool print_outputs = false;   // print_outputs on screen
    bool gif_style = false;       // save pictures from 1 to n and then again from n-1 to 1 to make nice .gif image later
    bool gt = false;              // ground truth given ?
    int MAX_RES = 100000;         // Max resolution for images. If reached, images are reduces to this resolution
    bool test = false;            // Test given method and saves in output dir a .csv file for plotting

    // variable du main
    bool first = true;







    // PROCESS ARGUMENTS
    // ===================================================
    char tmp;
    /*if the program is ran without options ,it will show the usgage and exit*/
    if(argc == 1){
        showhelpinfo(string(argv[0]));
        return EXIT_FAILURE;
    }
    // get options
    while ((tmp = getopt(argc, argv, "hvsgpd:r:o:m:e:t")) != -1){
      
      switch (tmp){
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
        /*option g for gif style format : output are saved from 0 to n then again from n-1 to 0 to make a double sens giff*/
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
        /*option e for error estimation with given */
        case 'e':
            ground_truth_path = string(optarg);
            gt = true;
            break;
        /*option t testing current method (output test_method.csv will appear in output_dir */
        case 't':
            test=true;
            break;
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
        cout << "Input directory option required (-d) !" << endl;
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

    // error evaluation : ground truth loading
    if(gt && ground_truth_path==""){
      cout << "Ground truth map required (after -e) !" << endl;
    }













    // CHARGEMENT DES IMAGES
    // ===================================================
    if(verbose){
        cout << endl;
        cout << "==================================================="<<endl;
        cout << "CHARGEMENT DES IMAGES" << endl;
        cout << "==================================================="<<endl;
    }
    // vector of images
    vector<Image<FVector<float,3>, 2 > > images;
    vector<Image<Color, 2 > > output_image_flow;
    vector<Image<float, 2 > > output_error_map;
    Image<FVector<float, 2>, 2 > ground_truth_map;

    //open directory and load images
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (directory.c_str())) != NULL) {
    /* print all the files and directories within directory */
        int count = 1;
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
                if(test && count++==2) break; // if test procedure, we only need the two first images
            }
            else{
                if(verbose) cout << "Ignoring uncompatible format file : " << name << endl;
            }
        }
        closedir (dir);
    }
    else {
    /* could not open directory */
        cout << string("Could not open directory ") + directory << endl;
        return EXIT_FAILURE;
    }

    // On vérifie le nombre d'images
    if(images.size()<2){
        cout << "Not enough images loaded, need at least 2" << endl;
        return EXIT_FAILURE;
    }



    // now load ground truth image
    if(gt){
      // verify that is is a .flo file
      if( ground_truth_path.find(".flo")==string::npos ){
        cout << "Error : ground truth map must be a .flo file !" << endl;
        return EXIT_FAILURE;
      }
      // load ground truth map
      try{
        ground_truth_map = flow_from_file(ground_truth_path);
        if(verbose) cout << "Ground truth correctly loaded" << endl;
      }
      catch(string e){
          cout << e;
        return EXIT_FAILURE;
      }

    }











    // PRETRAITEMENT DES IMAGES
    // ===================================================
    if(verbose){
        cout << endl;
        cout << "==================================================="<<endl;
        cout << "PRETRAITEMENT DES IMAGES" << endl;
        cout << "==================================================="<<endl;
    }

    int w = images[0].width(), h = images[0].height();
    bool has_been_pretraited = false;

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
        has_been_pretraited = true;
      }

    }

    // also verify that ground truth map
    if(gt){
        if(w != ground_truth_map.width() || h != ground_truth_map.height()){
            // Erreur si images de taille différentes
            throw string("Erreur : Ground truth et images de tailles différentes");
        }
        if(ground_truth_map.width() * ground_truth_map.height() > MAX_RES){
            // Rescale si on depasse la resolution maximale
            double fact = 1.0*MAX_RES / (ground_truth_map.width() * ground_truth_map.height());
            if(verbose) cout << "Reducing image ground truth map" << endl;
            ground_truth_map = reduce(ground_truth_map, 1/fact);
            has_been_pretraited = true;
      }
    }


    if(!has_been_pretraited && verbose) cout << "None" << endl;


    // affichage de la groud_truth
    if(gt){
      Image<Color, 2 > ground_truth_image = make_flow_visible_hsv(ground_truth_map);
      
      if(first && print_outputs){
        openWindow(ground_truth_image.width(), ground_truth_image.height());
        first=false;
      }
      if(print_outputs){
        display(ground_truth_image);
      }
    }








    // TRAITEMENT DES IMAGES
    // ===================================================
    if(verbose){
        cout << endl;
        cout << "==================================================="<<endl;
        cout << "TRAITEMENT DES IMAGES ..." << endl;
        cout << "==================================================="<<endl;
    }

    vector<pair<string, vector<string> > > test_results; // utilisé pour les tests : paires nom de fichier .csv et lignes à écrire dans ce fichier .csv qui pourra ensuite être tracé (Excel, Gnuplot, etc..)

    /*
        procedure normale : on calcule le flot optique pour les couples d'images prises 2 à 2
    */
    if(!test){
      // boucle sur les images
      for(int i=0; i<images.size()-1; i++){

        // Calcul du flow optique
        // ==============================
        if(verbose) cout << "Processing image n°" << i << endl;

        Image<FVector<float,2> ,2 > optical_flow;
      	if (method_args[0] == "LK"){
      		int taille_fenetre = atoi(method_args[1].c_str());    // Ne marche pas pour le moment car la taille de la fenêtre doit être fixé avant compilation, ici 7 (cf FVector)
      		optical_flow = flow_Lucas_Kanade(images[i], images[i + 1], taille_fenetre);
      	}
      	else if (method_args[0] == "HS"){
      		float smoothness = atof(method_args[1].c_str());
      		float stop = atof(method_args[2].c_str());
      		int max_iter = atoi(method_args[3].c_str());
      		optical_flow = flow_Horn_Schunk(images[i], images[i + 1], smoothness, stop, max_iter);
      	}
      	else if (method_args[0] == "HSL1"){
      		int max_iter = atoi(method_args[1].c_str());
      		optical_flow = flow_Horn_Schunk_HuberL1(images[i], images[i + 1], max_iter);
      	}
        // ==============================


        // Visualisation
        // ==============================
        Image<Color, 2 > optical_flow_image = make_flow_visible_hsv(optical_flow);
        output_image_flow.push_back(optical_flow_image);
        
        if(first && print_outputs){
          openWindow(optical_flow_image.width(), optical_flow_image.height()); 
          first=false;
        }
        if(print_outputs){
          display(optical_flow_image);
          anyClick();
        }
        // ==============================


        // Calcul de la map d'erreurs
        // ==============================
        if(gt){
            Image<float, 2> err_map = error_map(optical_flow, ground_truth_map );
            //err_map/=100.0; // rescale for not having overflow
            cout << "error is : " << sum(err_map) << endl;
        }
        // ==============================

      }
    }

    /*
        procédure de test. On teste une méthode pour plusieurs paramètres sur les deux premières images,
        et on renvoie un fichier .csv pour visualiser une courbe d'erreur en fonction des paramètres
    */
    else{
        Image<FVector<float,2> ,2 > optical_flow; // flow optique à calculer

        if(verbose) cout << "Testing method " << method_args[0] << endl;



        // Tester Luka et Kanade
        // ==============================
        vector<string> lines;
        int max_taille_fenetre = 19;
        // boucler sur la taille des fenêtres
        for(int taille_fenetre=3; taille_fenetre<=max_taille_fenetre; taille_fenetre+=2){

            if(verbose) cout << "test num " << taille_fenetre << " out of " <<  max_taille_fenetre << endl;
            optical_flow = flow_Lucas_Kanade(images[0], images[1], taille_fenetre);
            
            // calculer la carte d'erreur
            Image<float, 2> err_map = error_map(optical_flow, ground_truth_map );
            float sum_error = sum(err_map);
            if(verbose) cout << "LK "<< taille_fenetre<< "  |   error is : " << sum_error << endl << endl;
            lines.push_back( std::to_string(taille_fenetre) + ";" + std::to_string(sum_error) );
        }
        // add lines to test_results to be written
        test_results.push_back( pair<string, vector<string> >(output_directory+"/"+"LK_test_taille_fenetre.csv", lines) );
        // ==============================



        // Tester Horn et Schunk itératif
        // ==============================
        vector<string> lines;
        float max_smoothness = 200.0;
        float max_stop = 1.0;
        int max_iter = 100;

        // boucle sur les valeurs de smoothness avec stop fixée
        for(float smoothness=10.0; smoothness<=max_smoothness; smoothness+=5.0){

            if(verbose) cout << "test num " << smoothness/10 << " out of " <<  max_smoothness/10 << endl;
            optical_flow = flow_Horn_Schunk(images[0], images[1], smoothness, 0.02, max_iter);
            
            // calculer la carte d'erreur
            Image<float, 2> err_map = error_map(optical_flow, ground_truth_map );
            float sum_error = sum(err_map);
            if(verbose) cout << "HS smoothness "<< smoothness << "  |   error is : " << sum_error << endl << endl;
            lines.push_back( std::to_string(smoothness) + ";" + std::to_string(sum_error) );
        }
        // add lines to test_results to be written later
        test_results.push_back( pair<string, vector<string> >(output_directory + "/" + "HS_test_smoothness.csv", lines) );

        // boucle sur les valeurs de stop avec smoothness fixée
        lines.clear();
        int count = 0;
        for(float stop=0.01; stop<=10.0; stop*=2.0){

            if(verbose) cout << "test num " << count << " out of " <<  std::to_string(13) << endl;
            optical_flow = flow_Horn_Schunk(images[0], images[1], 10.0, stop, max_iter);
            // calculer la carte d'erreur
            Image<float, 2> err_map = error_map(optical_flow, ground_truth_map );
            float sum_error = sum(err_map);
            if(verbose) cout << "HS stop "<< stop << "  |   error is : " << sum_error << endl << endl;
            lines.push_back( std::to_string(stop) + ";" + std::to_string(sum_error) );
        }
        // add lines to test_results to be written later
        test_results.push_back( pair<string, vector<string> >(output_directory+"/"+"HS_test_stop.csv", lines) );
        // ==============================


        // Pas encore de tests pour HuberL1 (non fonctionnel)
        // ==============================
    }













    // SAUVEGARDE DES RESULTATS
    // ===================================================
    if(save_outputs){
      cout << endl;
      cout << "==================================================="<<endl;
      cout << "SAUVEGARDE DES RESULTATS" << endl;
      cout << "==================================================="<<endl;
    }

    if (save_outputs || test){

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
    }

    if( save_outputs ){

      if(verbose) cout << "Saving images to " + output_directory << endl;
      int n_dec = int(log10(output_image_flow.size()) + 1); // number of decimal to use for number in name
      char* num = new char[n_dec]; // tableau contenant le numéro  de l'image en chaine de caractères
      for(int i=0; i<output_image_flow.size(); i++){
        cout << "saving " << "/output_" << to_string(i) << ".jpg" << endl;
          
        sprintf(num, (string("%0"+to_string(n_dec)+"d")).c_str(), i);
        save(output_image_flow[i], output_directory + "/output_" + num +".jpg");
      }
      if(gif_style){ // save als in reverse for making a nice .gif image
        int j=output_image_flow.size();
        for(int i=output_image_flow.size()-1; i>=0; i--){ 
          cout << "saving " << "output_" << to_string(j) << ".png" << endl;
          save(output_image_flow[i], output_directory + "/output_" + to_string(j)+".png");
          j++;
        }
      }
    }

    // If test procedure, save csv files
    if ( test ){

        for(vector<pair<string, vector<string> > >::iterator it=test_results.begin();
            it!=test_results.end();
            it++){
            write_lines_to_csv( it->first , it->second );
        if(verbose) cout << "Wrote file " << it->first << endl;
        }

    }


return EXIT_SUCCESS;
}


/*
    fonction that shows the help information
*/
void showhelpinfo(string s)
{
    cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
    cout<<"option:  "<<"-h  show help information"<<endl;
    cout<<"         "<<"-d directory"<<endl;
    cout<<"         "<<"-o output directory"<<endl;
    cout<<"         "<<"-s save"<<endl;
    cout<<"         "<<"-t test given method"<<endl;
    cout<<"         "<<"-p print results"<<endl;
    cout<<"         "<<"-m method args \"LK taille_fenetre\"  or\n           \"HS smoothness stop max_iter\"  or\n        \"HSL1 max_iter\" "<<endl;
    cout<<"         "<<"-v verbose"<<endl;
    cout<<"         "<<"-e error estimation with given ground-truth"<<endl;
    cout<<"         "<<"-r max resolution of ouptputs"<<endl;
    cout<<"         "<<"-g gif_style"<<endl;
    cout<<"example: "<<s<<" -d directory -s1"<<endl;
}


/*
    writes lines to given file_path
*/
void write_lines_to_csv(const string& file_path, vector<string>& lines){
    ofstream file; // out file stream
    file.open(file_path);
    for(vector<string>::iterator it=lines.begin(); it!=lines.end(); it++){
        file << *it << endl;
    }
    file.close();
}
