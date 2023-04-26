#include <cmath>
#include <chrono>
#include <algorithm>
#include <iostream>
#include "state.h"

/*!
 * \brief Constructor of State
 *
 * Initializes the state of the system: particles randomly placed in a 2d box.
 */
State::State(const double _len_x, const double _len_y, const long _n_parts,
		      const double _pot_strength, const double _dt) :
	len_x(_len_x), len_y(_len_y), n_parts(_n_parts),
	pot_strength(_pot_strength), dt(_dt)
#ifdef USE_MKL
#else
	, rng(std::chrono::system_clock::now().time_since_epoch().count())
#endif
{
	positions[0].resize(n_parts);
	positions[1].resize(n_parts);
	forces[0].assign(n_parts, 0);
	forces[1].assign(n_parts, 0);

#ifdef USE_MKL
	aux_x.resize(n_parts);
	aux_y.resize(n_parts);
	vslNewStream(&stream, VSL_BRNG_SFMT19937,
		std::chrono::system_clock::now().time_since_epoch().count());
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
			     positions[0].data(), 0, len_x);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
			     positions[1].data(), 0, len_y);
#else
    std::uniform_real_distribution<double> rnd(0, 1.);

	for (long i = 0 ; i < n_parts ; ++i) {
		positions[0][i] = len_x * rnd(rng);
		positions[1][i] = len_y * rnd(rng);
		forces[0][i] = 0;
		forces[1][i] = 0;
	}
#endif
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolve() {
	calcInternalForces();

	for (long i = 0 ; i < n_parts ; ++i) {
		// Internal forces
		positions[0][i] += dt * forces[0][i];
		positions[1][i] += dt * forces[1][i];
	}

	enforcePBC();
}

void State::dump() const {
	for (long i = 0 ; i < n_parts ; ++i) {
		std::cout << positions[0][i] << " " << positions[1][i] << "\n";
	}
}

/*
 * \brief Compute the forces between the particles.
 */
void State::calcInternalForces() {
    for (long i = 0 ; i < n_parts ; ++i) {
		forces[0][i] = 0;
		forces[1][i] = 0;
    }

	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; i < j ; ++i) {
			calcInternalForceIJ_WCA(i, j);
		}
	}
}

//! Compute internal force between particles i and j (WCA potential)
void State::calcInternalForceIJ_WCA(const long i, const long j) {
	//std::cout << i << " " << j << "\n";

	double dx = positions[0][i] - positions[0][j];
	double dy = positions[1][i] - positions[1][j];
	// We want the periodized interval to be centered in 0
	pbcSym(dx, len_x);
	pbcSym(dy, len_y);
	double dr2 = dx * dx + dy * dy;

	if(dr2 * (TWOONESIXTH - dr2) > 0.) {
		double u = pot_strength;
		u *= (48. * pow(dr2, -7.) - 24.*pow(dr2, -4.)); 
		double fx = u * dx;
		double fy = u * dy;

		forces[0][i] += fx;
		forces[0][j] -= fx;
		forces[1][i] += fy;
		forces[1][j] -= fy;
	}
}

/* 
 * \brief Enforce periodic boundary conditions
 */
void State::enforcePBC() {
#ifdef USE_MKL
	pbcMKL(positions[0], len_x, aux_x, n_parts);
	pbcMKL(positions[1], len_y, aux_y, n_parts);
#else
	for (long i = 0 ; i < n_parts ; ++i) {
		pbc(positions[0][i], len_x);
		pbc(positions[1][i], len_y);
	}
#endif
}

#ifdef USE_MKL
void pbcMKL(std::vector<double> &v, const double L, std::vector<double> &aux,
	        const long N) {
	cblas_daxpby(N, 1.0 / L, v.data(), 1, 0.0, aux.data(), 1);
	vdFloor(N, aux.data(), aux.data());
	cblas_daxpy(N, -L, aux.data(), 1, v.data(), 1);
}	

void pbcSymMKL(std::vector<double> &v, const double L,
		       std::vector<double> &aux, const long N) {
	cblas_daxpby(N, 1.0 / L, v.data(), 1, 0.0, aux.data(), 1);
	vdRound(N, aux.data(), aux.data());
	cblas_daxpy(N, -L, aux.data(), 1, v.data(), 1);
}	
#endif
