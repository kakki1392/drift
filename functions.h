#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include "system.h"
#include <string>

using namespace std;

double potential(double &x, const double &alpha);

double delx_potential(double &x, const double &alpha);

double flashing(double &t, const double &tau, const bool &on);

double findTimeStep(const double &D, const double &alpha, double &k);

gsl_histogram * unitGaussHistogram(int N, int bin, double xmin, double xmax, const gsl_rng_type *T, int seed);

double nextEulerPos(const System &system, double &x, double &t, double &dt, double &square, gsl_rng * rng);

void positionTimeSeries(const System &system, std::string &filename, double &simTime_SI, double &k, gsl_rng * rng);

void positionTimeEndpoint(const System &system, double &simTime_SI, double &k, gsl_rng * rng, double &x, double &t);

double driftVelocity(const System &system, double &simTime_SI, double &k, gsl_rng * rng);

void driftVelocityStatistics(const System &system, double &simTime_SI, double &k, gsl_rng * rng, double &v_mean, double &v_std, size_t &iterations);

void driftVelocityStatisticsSeries(const System &system, double &simTime_SI, double &k, gsl_rng * rng, double &tau_min, double &tau_max, size_t &num_tau, size_t &iterations, string &filename);

void particleDensityMovement(const System &system, double &simTime_SI, double &k, size_t &num_simTime, size_t &num_particles, size_t &bins, gsl_rng * rng, std::string &filename);

#endif
