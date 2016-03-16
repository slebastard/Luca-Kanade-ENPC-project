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
    string directory("");         // Dossier d'inputs
    string output_directory("");  // Dossier d'outputs
    string ground_truth_path(""); // Chemin vers la ground_truth_map
    vector<string> method_args;   // arguments pour l'option -m (method) à écrire entre ""
    bool verbose = false;         // verbose
    bool save_outputs = false;    // sauvegarder les images
    bool print_outputs = false;   // afficher les images
    bool gif_style = false;       // Sauvegarder les images aussi en ordre inverse pour faire de beaux .gif
    bool gt = false;              // Ground truth donnée ou non
    int MAX_RES = 100000;         // Résolution max pour les images, si atteinte, les images sont réduites en prétraitement
    bool test = false;            // Permet de tester les méthodes et écrit un fichier .csv pour tracer des courbes

    // variable du main
    bool first = true;





    // TRAITEMENT DES OPTIONS PASSEES SUR LA LIGNE DE COMMANDE
    // ===================================================
    char tmp;
    // Si pas d'option, afficher l'aide et uitter
    if(argc == 1){
        showhelpinfo(string(argv[0]));
        return EXIT_FAILURE;
    }
    
    while ((tmp = getopt(argc, argv, "hvsgpd:r:o:m:e:t")) != -1){
      
      switch (tmp){
        /* Option h affiche l'aide */
        case 'h':
            showhelpinfo(string(argv[0]));
            break;
        /* Option v pour verbose*/
        case 'v':
            verbose = true;
            break;
        /* Option s pour sauvegarder */
        case 's':
            save_outputs = true;
            break;
        /* Option r pour la resolution max */
        case 'r':
            MAX_RES = atoi(optarg);
            break;
        /* Option o pour le dossier d'outputs */
        case 'o':
            output_directory = string(optarg);
            break;
        /* Option g pour enregistre aussi les images en sens inverse pour faire de beaux .gif */
        case 'g':
            gif_style = true;
            break;
        /* Option p pour afficher les images */
        case 'p':
            print_outputs = true;
            break;
        /* Option d pour le dossier d'inputs */
        case 'd':
            directory = string(optarg);
            break;
        /* Option e pour l'emplacement de la ground_truth_map (.flo) */
        case 'e':
            ground_truth_path = string(optarg);
            gt = true;
            break;
        /* Option t pour lancer les tests */
        case 't':
            test=true;
            break;
        break;
        /* Option m pour spécifier une méthode de calcul de flux optique */
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

    // On vérifie les options passées

    // directory
    if(directory==""){
        cout << "Input directory option required (-d) !" << endl;
        return EXIT_FAILURE;
    }
    if(output_directory==""){
        output_directory = directory + "/outputs";
    }

    // save_outputs
    if(!save_outputs) print_outputs = true;

    // method
    if(!test && method_args.size()==0){
        cout << "Method option required !" << endl;
        return EXIT_FAILURE;
    }

    if (!test && method_args[0] != "LK" && method_args[0] != "HS" && method_args[0] != "HSL1"){
        cout << "Unknown method : " << method_args[0] << endl;
        showhelpinfo();
        return EXIT_FAILURE;
    }
    if(!test && method_args[0]=="LK"){
      if(method_args.size()!=2){
        cout << "Incorrect args number for LK method ! " << endl;
        showhelpinfo();
        return EXIT_FAILURE;
      }
    }
    else if(!test && method_args[0]=="HS"){
      if(method_args.size()!=4){
        cout << "Incorrect args number for HS method ! " << endl;
        showhelpinfo();
        return EXIT_FAILURE;
      }
    }
    else if (!test && method_args[0] == "HSL1"){
    	if (method_args.size() != 2){
    		cout << "Incorrect args number for HSL1 method ! " << endl;
        showhelpinfo();
    		return EXIT_FAILURE;
    	}
    }

    // Chargement de la ground_truth
    if(gt && ground_truth_path==""){
      cout << "Ground truth map required (after -e) !" << endl;
    }













    // CHARGEMENT DES IMAGES
    // ===================================================
    if(verbose){
        cout << endl;
        cout << "==================================================="<<endl;
        cout << "LOADING IMAGES" << endl;
        cout << "==================================================="<<endl;
    }
    // Vecteur d'images
    vector<Image<FVector<double,3>, 2 > > images;
    vector<Image<Color, 2 > > output_image_flow;
    vector<Image<double, 2 > > output_error_map;
    Image<FVector<double, 2>, 2 > ground_truth_map;

    // Ouvre le dossier d'inputs et charge les images
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (directory.c_str())) != NULL) {
    // Affiche les fichiers trouvés
        int count = 1;
        while ((ent = readdir (dir)) != NULL) {
            string name(ent->d_name);
            // Recherche des formats png jpg et tiff
            if(name.find(".png")!=string::npos || name.find(".jpg")!=string::npos || name.find(".tiff")!=string::npos ){
                Image<Color> I;
                if(!load(I, directory + "/" + name)){
                    cout << string("Couldn't load image : ") + name << endl;
                    return EXIT_FAILURE;
                }
                Image<FVector<double,3> > Iv(I);
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
        cout << string("Could not open directory ") + directory << endl;
        return EXIT_FAILURE;
    }

    // On vérifie le nombre d'images
    if(images.size()<2){
        cout << "Not enough images loaded, need at least 2" << endl;
        return EXIT_FAILURE;
    }



    // On charge la ground_truth si demandé
    if(gt){
      // Vérifier qu'il s'agit d'un .flo
      if( ground_truth_path.find(".flo")==string::npos ){
        cout << "Error : ground truth map must be a .flo file !" << endl;
        return EXIT_FAILURE;
      }
      // charger la ground_truth
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
        if(verbose) cout << "Reduction de l'image n°" << i << endl;
        images[i] = reduce(images[i], 1/fact);
        has_been_pretraited = true;
      }

    }

    // Verifier aussi la ground truth map
    if(gt){
        if(w != ground_truth_map.width() || h != ground_truth_map.height()){
            // Erreur si images de taille différentes
            throw string("Erreur : Ground truth et images de tailles différentes");
        }
        if(ground_truth_map.width() * ground_truth_map.height() > MAX_RES){
            // Rescale si on depasse la resolution maximale
            double fact = 1.0*MAX_RES / (ground_truth_map.width() * ground_truth_map.height());
            if(verbose) cout << "Réduction de la ground truth map" << endl;
            ground_truth_map = reduce(ground_truth_map, 1/fact);
            has_been_pretraited = true;
      }
    }


    if(!has_been_pretraited && verbose) cout << "None" << endl;


    // Affichage de la groud_truth
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
        if(verbose) cout << "Traitement de l'image n°" << i << endl;

        Image<FVector<double,2> ,2 > optical_flow;
      	if (method_args[0] == "LK"){
      		int taille_fenetre = atoi(method_args[1].c_str());
      		optical_flow = flow_Lucas_Kanade(images[i], images[i + 1], taille_fenetre);
      	}
      	else if (method_args[0] == "HS"){
      		double smoothness = atof(method_args[1].c_str());
      		double stop = atof(method_args[2].c_str());
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
            Image<double, 2> err_map = error_map(optical_flow, ground_truth_map );
            //err_map/=100.0; // rescale for not having overflow
            cout << "Difference avec la ground_truth : " << sum(err_map) << endl;
        }
        // ==============================

      }
    }

    /*
        procédure de test. On teste une méthode pour plusieurs paramètres sur les deux premières images,
        et on renvoie un fichier .csv pour visualiser une courbe d'erreur en fonction des paramètres
    */
    else{
        Image<FVector<double,2> ,2 > optical_flow; // flow optique à calculer

        if(verbose) cout << "Méthode testée : " << method_args[0] << endl;



        // Tester Luka et Kanade
        // ==============================
        vector<string> lines;
        int max_taille_fenetre = 19;
        // boucler sur la taille des fenêtres
        for(int taille_fenetre=3; taille_fenetre<=max_taille_fenetre; taille_fenetre+=2){

            if(verbose) cout << "test numéro " << taille_fenetre << " sur " <<  max_taille_fenetre << endl;
            optical_flow = flow_Lucas_Kanade(images[0], images[1], taille_fenetre);
            
            // calculer la carte d'erreur
            Image<double, 2> err_map = error_map(optical_flow, ground_truth_map );
            double sum_error = sum(err_map);
            if(verbose) cout << "LK "<< taille_fenetre<< "  |   erreur : " << sum_error << endl << endl;
            lines.push_back( std::to_string(taille_fenetre) + ";" + std::to_string(sum_error) );
        }
        // Ajouter les lignes à écrire à test_results
        test_results.push_back( pair<string, vector<string> >(output_directory+"/"+"LK_test_taille_fenetre.csv", lines) );
        // ==============================


        // Tester Horn et Schunk itératif
        // ==============================
        double max_smoothness = 10.0f;
        double max_stop = 0.2f;
        int max_iter = 100;

        // Boucle sur les valeurs de smoothness avec stop fixé
        lines.clear();
        int count = 0;
        for(double smoothness=0.5; smoothness<=max_smoothness; smoothness+=0.5){
            for(double stop=0.001; stop<=max_stop; stop*=2.0){

                if(verbose) cout << "test numéro " << count++ << " sur " <<  std::to_string((int)13 * max_smoothness/10) << endl;
                optical_flow = flow_Horn_Schunk(images[0], images[1], smoothness, stop, max_iter);
                // calculer la carte d'erreur
                Image<double, 2> err_map = error_map(optical_flow, ground_truth_map );
                double sum_error = sum(err_map);
                if(verbose) cout << "HS smoothness " <<  smoothness << " & stop "<< stop << "  |   erreur : " << sum_error << endl << endl;
                lines.push_back( std::to_string(smoothness) + ";" + std::to_string(stop) + ";" + std::to_string(sum_error) );
            }
        }
        // Ajouter les lignes à écrire à test_results
        test_results.push_back( pair<string, vector<string> >(output_directory+"/"+"HS_test.csv", lines) );
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
      int n_dec = int(log10(output_image_flow.size()) + 1);
      char* num = new char[n_dec]; // tableau contenant le numéro de l'image en chaine de caractères
      for(int i=0; i<output_image_flow.size(); i++){
        cout << "saving " << "/output_" << to_string(i) << ".jpg" << endl;
          
        sprintf(num, (string("%0"+to_string(n_dec)+"d")).c_str(), i);
        save(output_image_flow[i], output_directory + "/output_" + num +".jpg");
      }
      if(gif_style){ // si option -g spécifiée
        int j=output_image_flow.size();
        for(int i=output_image_flow.size()-1; i>=0; i--){ 
          cout << "saving " << "output_" << to_string(j) << ".png" << endl;
          save(output_image_flow[i], output_directory + "/output_" + to_string(j)+".png");
          j++;
        }
      }
    }

    // Si procédure de test, sauvegarder les résultats en .csv
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
    Fonction d'aide
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
    Ecrit des lignes dans un fichier csv
*/
void write_lines_to_csv(const string& file_path, vector<string>& lines){
    ofstream file; // out file stream
    file.open(file_path);
    for(vector<string>::iterator it=lines.begin(); it!=lines.end(); it++){
        file << *it << endl;
    }
    file.close();
}
