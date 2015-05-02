#ifndef SYSTEM_H
#define SYSTEM_H
#include <string>

const double eV = 1.602176565e-19;
const double PI = 3.14159265359;

class System {
	public:
		System();
		System(std::string filename);
		System(const System &system);
		~System();
		void writeSystemParameters(std::string filename);
		void printSystemParameters();
		void initializeSystem(std::string filename);
		void copy(const System &system);	

		double getTau();
		void setTau(double tau_new);
		void updateTau(double d_tau);

		double r_SI;
		double L_SI;
		double alpha;
		double kbT_eV;
		double deltaU_eV;
		double tau_SI;
		double eta_SI;
		
		bool isFlashing;

		double omega;
		double gamma;
		double D;
		double tau;

		void makeDimensionless();
};





#endif


