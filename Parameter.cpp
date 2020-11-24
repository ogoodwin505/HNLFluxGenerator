
#include "Parameter.h"


#include <iostream>
#include <iomanip>
#include <string>






#include <fstream>



#include "TChain.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TRotation.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"

#include <cmath>
#include <vector>
#include <map>
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

//___________________________________________________________________________
int split_string(std::string &hold, std::string &key, std::string &val) {
  if (hold.length() == 0 || hold[0] == '#')
    return 0;

  std::size_t space_pos;
  if ((space_pos = hold.find(" ")) == std::string::npos)
    return -1;

  key = hold.substr(0, space_pos);
  val = hold.substr(space_pos + 1, hold.find(" ", space_pos + 1));
  return 1;
}


Parameter::Parameter(const std::string &filename) {

  std::ifstream instream(filename);
  if (instream.is_open()) {
    std::string hold, key, val;
    int line_num = 0;
    int error_state = 0;
    integrity=0;
    //Must be used if you are using different parsing modes (like materials!)
    //int last_pos = instream.tellg();
    while (getline(instream, hold)) {
      line_num++;
      if ((error_state = split_string(hold, key, val)) == 0);
      else if (error_state == -1) {
        std::cerr << line_num << ": Improperly formatted input: ";
        std::cerr << hold << "\n";
      }
      else if (error_state != 1) {
        std::cerr << line_num << ": Unknown error code: " << error_state << "\n";
      }
      else {
        keymap[key] = val;
      }
      //last_pos=instream.tellg();
    }

    instream.close();
  }
  else{
    std::cerr << "Could not open " << filename << std::endl;
  }
  Set_Double("det_posx", det_posx, keymap);
  Set_Double("det_posy", det_posy, keymap);
  Set_Double("det_posz", det_posz, keymap);


  Set_Double("TPC_hix",TPC_hix, keymap);
  Set_Double("TPC_lox",TPC_lox, keymap);
  Set_Double("TPC_hiy",TPC_hiy, keymap);
  Set_Double("TPC_loy",TPC_loy, keymap);
  Set_Double("TPC_hiz",TPC_hiz, keymap);
  Set_Double("TPC_loz",TPC_loz, keymap);


  Set_String("outputloc", outputloc, keymap); 
}//end of Parameter


void Parameter::Set_Double(const std::string &key, double &var, std::map<std::string, std::string> &keymap){
    try{
        if(keymap.count(key)==1){
            var = stod(keymap[key]);
        }
        else{
            std::cerr << key << " not set." <<std::endl;
            integrity = -1;
        }
    }
    catch(std::exception& e){
            std::cerr << "Invalid model parameter: " << key << std::endl;
            integrity = -1;
    }
}


void Parameter::Set_String(const std::string &key, std::string &var, std::map<std::string, std::string> &keymap){
  if(keymap.count(key)==1)
    var = keymap[key];
  else
    var = "./"; //default if set
}



//___________________________________________________________________________

TVector3 Parameter::DetPos_vec() {
  TVector3 DetLoc (det_posx, det_posy, det_posz); //cm
  return DetLoc;
}

TVector3 Parameter::DetUpperBound() {
  TVector3 UpperBound (TPC_hix, TPC_hiy, TPC_hiz); //cm
  return UpperBound;
}


TVector3 Parameter::DetLowerBound() {
  TVector3 LowerBound (TPC_lox, TPC_loy, TPC_loz); //cm
  return LowerBound;
}

