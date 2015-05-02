#include "system.h"
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace std;

void System::makeDimensionless(){
	D = kbT_eV/deltaU_eV;
	gamma = 6.0*PI*eta_SI*r_SI;
	omega = deltaU_eV*eV/(gamma*L_SI*L_SI);
	tau = omega*tau_SI;
}

System::System(){
	r_SI = 12.0e-9;
	L_SI = 20.0e-6;
	alpha = 0.2;
	kbT_eV = 26.0e-3;
	deltaU_eV = 80.0;
	isFlashing = false;
	tau_SI = 1.0;
	eta_SI = 1.0e-3;
	makeDimensionless();
}

System::System(std::string filename){
	initializeSystem(filename);
	makeDimensionless();
}

void System::copy(const System &system){
	r_SI = system.r_SI;
	L_SI = system.L_SI;
	alpha = system.alpha;
	kbT_eV = system.kbT_eV;
	deltaU_eV = system.deltaU_eV;
	tau_SI = system.tau_SI;
	eta_SI = system.eta_SI;
	isFlashing = system.isFlashing;
	
	makeDimensionless();
}

System::System(const System &system){
	copy(system);
}

System::~System(){};



void System::initializeSystem(std::string filename){
	ifstream input;
	input.open(filename.c_str());
	string line, name, number;
	size_t pos1, pos2, pos3;
	while(getline(input,line)){
		if(line[0] == '#' || line[0] == ' '){
			continue;
		}

		pos1 = line.find(' ');
		pos2 = line.find(' ', pos1+1);
		pos3 = line.find(' ', pos2+1);
		name = line.substr(0,pos1);
		number = line.substr(pos2+1, pos3-pos2-1);
		if(name == "radius"){
			r_SI = atof(number.c_str());
		}else if(name == "length"){
			L_SI = atof(number.c_str()); 
		}else if(name == "alpha"){
			alpha = atof(number.c_str()); 
		}else if(name == "KbT"){
			kbT_eV = atof(number.c_str()); 
		}else if(name == "deltaU"){
			deltaU_eV = atof(number.c_str()); 
		}else if(name == "eta"){
			eta_SI = atof(number.c_str()); 
		}else if(name == "flash"){
			isFlashing = atoi(number.c_str()); 
		}else if(name == "tau"){
			tau_SI = atof(number.c_str()); 
		}else{
			continue;
		}
	}
	input.close();
}

void System::writeSystemParameters(string filename){
	ofstream output;
	output.open(filename.c_str());
	output << "#These are the physical parameters of the system in consideration\n\n";
	output << "#Radius of particle: \n";
	output << "radius = " << r_SI << " m\n\n";
	output << "#Size of potential: \n";
	output << "length = " << L_SI << " m\n\n";
	output << "#Asymmetri factor: \n";
	output << "alpha = " << alpha << " \n\n";
	output << "#Thermal energy: \n";
	output << "KbT = " << kbT_eV << " eV\n\n";
	output << "#Strength of potential: \n";
	output << "deltaU = " << deltaU_eV << " eV\n\n";
	output << "#Strength of friction: \n";
	output << "eta = " << eta_SI << " eV\n\n";
	output << "#Is the potential flashing?: \n";
	output << "flash = " << isFlashing << " \n\n";
	output << "#Flashing period of potential: \n";
	output << "tau = " << tau_SI << " s\n"; 

	output.close();
}

void System::printSystemParameters(){
	cout << "#These are the physical parameters of the system in consideration\n \n";
	cout << "#Radius of particle: \n";
	cout << "radius = " << r_SI << " m\n\n";
	cout << "#Size of potential: \n";
	cout << "length = " << L_SI << " m\n\n";
	cout << "#Asymmetri factor: \n";
	cout << "alpha = " << alpha << " \n\n";
	cout << "#Thermal energy: \n";
	cout << "KbT = " << kbT_eV << " eV\n\n";
	cout << "#Strength of potential: \n";
	cout << "deltaU = " << deltaU_eV << " eV\n\n";
	cout << "#Strength of friction: \n";
	cout << "eta = " << eta_SI << " sPa \n\n";
	cout << "#Is the potential flashing?: \n";
	cout << "flash = " << isFlashing << " \n\n";
	cout << "#Flashing period of potential: \n";
	cout << "tau = " << tau_SI << " s\n"; 

}

double System::getTau(){ return tau_SI;}

void System::setTau(double tau_new){
	tau_SI = tau_new;
	tau = omega*tau_SI;
}

void System::updateTau(double d_tau){
	tau_SI = tau_SI + d_tau;
	tau = omega*tau_SI;
}







