#ifndef Parameter_h
#define Parameter_h

// #include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
// #include <list>
// #include <memory>

#include "TVector3.h"

class Parameter {
public:
	Parameter(const std::string &filename);
	~Parameter() {}


	TVector3 DetPos_vec();
	TVector3 DetUpperBound();
	TVector3 DetLowerBound();

	int Integrity(){return integrity;}

	double DetPos_x() {return det_posx;}
	std::string Output_loc() {return outputloc;}

private:
	std::map<std::string, std::string> keymap;


	std::string outputloc;
	//Detector Position in Beam Coords [cm]
	double det_posx;
	double det_posy;
	double det_posz;


	//Cuboid Detector Bounds [cm]
	double TPC_hix;
	double TPC_lox;
	double TPC_hiy;
	double TPC_loy;
	double TPC_hiz;
	double TPC_loz;



	void Set_Double(const std::string &, double &, std::map<std::string, std::string> &keymap);
	void Set_String(const std::string &, std::string &, std::map<std::string, std::string> &keymap);
	int integrity;
};

#endif
