#include "functions.h"
#include "gnuplotting.h"
#include <armadillo>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_randist.h>
#include <iostream>

using namespace std;
using namespace arma;

double potential(double &x, const double alpha){
	//This function returns a periodic asymmetrical sawtooth potential
	//with period 1. Asymmetri coefficient is alpha.

	double intpart, decimalpart; //separate x into integer+decimal
	decimalpart = modf(x, &intpart);
	if (decimalpart < 0.0){
		decimalpart = 1.0 + decimalpart; //Checks if we have negative x-values
	}

	if (decimalpart < alpha){
		return decimalpart/alpha;
	}else{
		return (1.0-decimalpart)/(1-alpha);
	}
}

double delx_potential(double &x, const double &alpha){
	//This function returns the derivative with respect to x of the function _double potential_
	double intpart, decimalpart;
	decimalpart = modf(x, &intpart);
	if (decimalpart < 0.0){
		decimalpart = 1.0 + decimalpart;
	}
	if (decimalpart < alpha){
		return 1/alpha;
	}else{
		return (-1.0)/(1-alpha);
	}
}

double flashing(double &t, const double &tau, const bool &on){
	//This function returns a value that is either 0 or 1.
	//It is 0 in 75% of the time. 
	//If on=true, the flashing is active. If not, the function always returns 1.0
	if(!on){
		return 1.0;
	}
	double t_eff = t/tau; //careful if large t! Need more robust implementation.
	double intpart, decimalpart;
	decimalpart = modf(t_eff, &intpart);
	if (decimalpart < (3.0/4.0)){
		return 0.0;
	}else{
		return 1.0;
	}
}

double findTimeStep(const double &D, const double &alpha, double &k){
	//Returns necessary timestep for a minimum of 1/k x-values in the region of length alpha*L ((1-alpha)*L), 
	//where L is period of potential.  0 < alpha < 1. 
	double alpha_temp = alpha;
	if (alpha > 0.5){
		alpha_temp = 1.0 - alpha;
	}
	return (16.0*D+k*alpha_temp - sqrt(16.0*16.0*D*D+32.0*D*k*alpha_temp));
}	


gsl_histogram * unitGaussHistogram(int N, int bin, double xmin, double xmax, const gsl_rng_type * T, int seed){
	gsl_rng *r = gsl_rng_alloc(T);
	gsl_rng_set(r,seed);
	gsl_histogram *h = gsl_histogram_alloc(bin);
	gsl_histogram_set_ranges_uniform(h,xmin,xmax);
	double x;
	for(int i=0; i<N; i++){
		x = gsl_ran_ugaussian(r);
		if( x<=xmin || x>=xmax){
			continue;
		}
		gsl_histogram_increment(h,x);
	}
	return h;
}

double nextEulerPos(const System &system, double &x, double &t, double &dt, double &square, gsl_rng * rng){
	//Returns the next position x at time t by stochastic Euler process. <System> contains all parameters of the physical system.
	//x, t, dt are all dimensionless.
	//gsl_rng*rng is a random number generator. 
	double xnew;
	xnew = x - delx_potential(x, system.alpha)*flashing(t, system.tau, system.isFlashing)*dt + square*gsl_ran_ugaussian(rng);
	return xnew;
}

void positionTimeSeries(const System &system, string &filename, double &simTime_SI, double &k, gsl_rng * rng){
	//Follows a particle, x=x(t), for a time simTime_SI seconds. Saves x(t) in dimensionless variables to
	//<filename>. The particle path is computed by the Euler process given in nextEulerPos().
	//Comments: -consider warning of possibly huge number of iterations.
	//- not saving all x(t) and t values.
       double dt = findTimeStep(system.D, system.alpha, k);	
       double simTime = system.omega*simTime_SI;
       size_t iterations = (size_t) (simTime/dt);

       ofstream output;
       output.open(filename.c_str());
       output << "#Position-Time series. Flashing period tau = " << system.tau_SI << " s" << endl;

       double x = 0.0;
       double t = 0.0;
       double square = sqrt(2.0*system.D*dt);
       output << t << " " << x << endl;
       for(size_t i = 0; i < iterations; i++){
	       x = nextEulerPos(system, x, t, dt, square, rng);
	       t = t + dt;
	       output << t/(system.omega) << " " << x << endl;
       }
       output.close();
}

void positionTimeEndpoint(const System &system, double &simTime_SI, double &k, gsl_rng * rng, double &x, double &t){
	//Follows a particle, x=x(t), for a time simTime_SI seconds. The last position and ending time is stored in x and t.
	//x and t are dimensionless.
	//Comments: -consider warning of possibly huge number of iterations.
       double dt = findTimeStep(system.D, system.alpha, k);	
       double simTime = system.omega*simTime_SI;
       size_t iterations = (size_t) (simTime/dt);
	
       x = 0.0;
       t = 0.0;
       double square = sqrt(2.0*system.D*dt);
       for(size_t i = 0; i < iterations; i++){
	       x = nextEulerPos(system, x, t, dt, square, rng);
	       t = t + dt;
       }
}

double driftVelocity(const System &system, double &simTime_SI, double &k, gsl_rng * rng){
	//Returns a drift velocity defined by (x_end - x_start)/(t_end - t_start). Here particle starts at the origin -> x_start=t_start=0.
	//A particle is followed in positionTimeEndpoint(), and the drift values are computed by the returned values of x and t of 
	//this function.
	double x, t;
	positionTimeEndpoint(system, simTime_SI, k, rng, x, t);
	return x/t;
}

void driftVelocityStatistics(const System &system, double &simTime_SI, double &k, gsl_rng * rng, double &v_mean, double &v_std, size_t &iterations){
	//Here a drift velocity is computed <iterations> times, and the mean and standard deviation are returned in v_mean and v_std.
	//The systems for each iteration are identical.
	double iterations_double = (double) iterations;
	v_mean = 0.0;
	double v_variance = 0.0;
	vector<double> v_vector (iterations);
	for(size_t i = 0; i < iterations; i++){
		double v_temp = driftVelocity(system, simTime_SI, k, rng);
		v_mean = v_mean + v_temp;
		v_vector[i] = v_temp;
	}
	v_mean = v_mean/iterations_double;
	for(size_t i = 0; i < iterations; i++){
		v_variance = v_variance + pow((v_vector[i] - v_mean), 2.0);
	}
	v_variance = v_variance/(iterations_double - 1.0);
	v_std = sqrt(v_variance);
}

void driftVelocityStatisticsSeries(const System &system, double &simTime_SI, double &k, gsl_rng * rng, double &tau_min, double &tau_max, size_t &num_tau, size_t &iterations, string &filename){
	//This function computes the mean and standard deviation of the drift velocity as a function of flashing time tau. There are num_tau tau-values
	//between tau_min and tau_max. For each tau-value a driftVelocityStatistics() with <iterations> run. A total of iterations*num_tau x(t) paths will have to be
	//computed, this may take a lot of time. tau-values and drift statistics are stored in <filename>.
	//Comments: add warning of long computation time. Can write an estimation of time required to console, with option to abort or proceed simulation.
	System system_iterator(system);
	system_iterator.setTau(tau_min);
	double d_tau = (tau_max - tau_min)/((double) (num_tau-1));

	ofstream output;
	output.open(filename.c_str());
	output << "#Drift velocities. Number of iterations per tau-value: " << iterations << ". Number of tau-values: " << num_tau << endl
		<< "#tau_min: " << tau_min << ". tau_max: " << tau_max << "." << endl << "# tau        v_mean               v_std" << endl;

	double v_mean = 0.0;
	double v_std = 0.0;

	for(size_t i = 0; i < num_tau; i++){
		cout << "tau number: " << i << endl;
		driftVelocityStatistics(system_iterator, simTime_SI, k, rng, v_mean, v_std, iterations);
		output << system_iterator.getTau() << " " << v_mean << " " << v_std << endl;
		system_iterator.updateTau(d_tau);
	}
	output.close();
}

void particleDensityMovement(const System &system, double &simTime_SI, double &k, size_t &num_simTime, size_t &num_particles, size_t &bins, gsl_rng * rng, string &filename){
	double simTime = system.omega*simTime_SI;
	double dt = findTimeStep(system.D, system.alpha, k);
	double square = sqrt(2.0*system.D*dt);
	size_t iterations = (size_t) (simTime/dt);	
	size_t extraction_point = (iterations-1)/num_simTime;

	mat pos_matrix (num_particles, num_simTime);
	for(size_t i = 0; i < num_particles; i++){
		double x = 0.0;
		double t = 0.0;
		size_t k = 0;
		for(size_t j = 0; j < iterations; j++){
			if( (j >= extraction_point) && (j % extraction_point == 0)){
				pos_matrix(i,k) = x;
				k = k + 1;
			}
			x = nextEulerPos(system, x, t, dt, square, rng);
			t = t + dt;
		}
		cout << "particle number: " << i+1 << " finished." << endl;
	}

	double time = 0.0;
	vector<double> time_points;
	for(size_t i = 0; i < iterations; i++){
		if((i >= extraction_point) && (i % extraction_point == 0)){
			time_points.push_back(time);
		}
		time = time + dt;
	}
	
	if(time_points.size() != (num_simTime)){
		cout << "ERROR IN GETTING ALL SIMULATION TIMES!!!!" << endl;
	}

	ofstream output;
	//output.open(filename.c_str());

	//add a delta function first
	/*
	for(size_t i = 0; i < bins; i++){
		output << "0.0 0.0 " << 1.0 << endl;
	}
	*/
	
	ofstream out;
	out.open("band.dat");
	for(size_t i = 0; i < (num_simTime); i++){
		for(size_t j = 0; j < num_particles; j++){
			out << time_points[i]/system.omega << " " << pos_matrix(j,i) << endl;
		}
		out << endl;
	}

	//make histograms
	Gnuplotting plot;
	cout << "starting binning and writing to file process" << endl;
	for(size_t i = 0; i < (num_simTime); i++){
		output.open(filename.c_str());
		vec bin_centers = linspace<vec>(pos_matrix.col(i).min(), pos_matrix.col(i).max(), bins);
		uvec histogram_int = hist(pos_matrix.col(i), bin_centers);
		vec histogram_double = conv_to< vec >::from(histogram_int);
		histogram_double = histogram_double / num_particles;
		for(size_t j = 0; j < bins; j++){
			output << time_points[i]/system.omega << " " << bin_centers(j) << " " << histogram_double(j) << endl;
		}
		//output.flush();
		output.close();
		cout << "Press 'e' for next plot: ";
		string str;
		cin >> str; 
		string command = "plot '" + filename + "' using 2:3 w impulses";
		plot.cmd(command);
	}
	cout << "finished writing" << endl;
}



