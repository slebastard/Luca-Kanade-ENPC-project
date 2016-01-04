/*
Projet MOPSI : Flux optique
Loic Cressot  &  Simon Lebastard
Encadrant : Pascal Monasse
janvier 2016
*/

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>


#include <Imagine/Images.h>
#include "optflow.hpp"


using namespace std;
using namespace Imagine;


void showhelpinfo(char *s);


int main (int argc,char *argv[])
{
    string directory("");
    bool verbose = false;

// PROCESS ARGUMENTS
    char tmp;
/*if the program is ran witout options ,it will show the usgage and exit*/
    if(argc == 1)
    {
        showhelpinfo(argv[0]);
        return EXIT_FAILURE;
    }
// get options
    while((tmp=getopt(argc,argv,"hvd:"))!=-1)
    {
        switch(tmp)
        {
  /*option h show the help infomation*/
          case 'h':
          showhelpinfo(argv[0]);
          break;
  /*option v for verbose*/
          case 'v':
          verbose = true;
          break;
  /*option d asks for directory*/
          case 'd':
          directory=string(optarg);
          break;
  // default shows help
          default:
          showhelpinfo(argv[0]);
          break;
      }
  }

// verify options
if(directory==""){
    cout << "Directory option required !" << endl;
    return EXIT_FAILURE;
}



// MAIN PROGRAMM

// vector of images
vector<Image<FVector<float,3> > > images;

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
            if(!load(I,directory + name)){
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

//checking images vector's size
if(images.size()<2){
    cout << "Not enough images loaded, need at least 2" << endl;
    return EXIT_FAILURE;
}


// At this point we should have a vector of Images<Color> of size 2 at least
find_flow_kanade(images[0], images[1]);



return 0;
}




/*funcion that show the help information*/
void showhelpinfo(char *s)
{
    cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
    cout<<"option:  "<<"-h  show help information"<<endl;
    cout<<"         "<<"-d directory"<<endl;
    cout<<"         "<<"-v verbose"<<endl;
    cout<<"example: "<<s<<" -d directory -s1"<<endl;
}
