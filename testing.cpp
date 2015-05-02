#include "functions.h"
#include "system.h"
#include "gnuplotting.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;


//HELLOOOOO



int main(){

System system;
double k = 0.1;
double sim = 5.0;
double period = 1.0;
system.tau_SI = period;
system.isFlashing = true;
//system.deltaU_eV = 0.000001;
system.makeDimensionless();

system.printSystemParameters();

const gsl_rng_type * T = gsl_rng_mt19937;
gsl_rng * rng = gsl_rng_alloc(T);
size_t seed = 10e16;
gsl_rng_set(rng, seed); 

string filename = "READ.dat";
//positionTimeSeries(system, filename, sim, k, rng);


double v_mean = 0.0;
double v_std = 0.0;

size_t iterations = 1000;
//driftVelocityStatistics(system,sim,k,rng,v_mean,v_std,iterations);

string file = "drift6.dat";

double tau_min = 0.01;
double tau_max = 1.5;
size_t num_tau = 40;
//driftVelocityStatisticsSeries(system,sim,k, rng, tau_min, tau_max, num_tau, iterations, file);

size_t num_sim = 3;
size_t num_part = 1000;
size_t bins = 10000;
string file2 = "density6.dat";
particleDensityMovement(system, sim, k, num_sim, num_part, bins, rng, file2);




	return 0;
}
