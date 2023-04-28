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
State::State(double _len_x, double _len_y, long _n_parts, double _a,
		     double _WCA_strength, double _dt, double _alpha_ew) :
	Lx(_len_x), Ly(_len_y), fac_x(1. / _len_x), fac_y(1. / _len_y),
	n_parts(_n_parts), sigma2(4 * _a * _a),
	WCA_strength(_WCA_strength), dt(_dt),
	ewald(_len_x, _len_y, _alpha_ew, _n_parts)
#ifdef USE_MKL
#else
	, rng(std::chrono::system_clock::now().time_since_epoch().count())
#endif
{
	positions[0].resize(n_parts);
	positions[1].resize(n_parts);
	forces[0].assign(n_parts, 0);
	forces[1].assign(n_parts, 0);
	dists[0].resize(n_parts * (n_parts - 1) / 2);
	dists[1].resize(n_parts * (n_parts - 1) / 2);

#ifdef USE_MKL
	aux_x.resize(n_parts);
	aux_y.resize(n_parts);
	vslNewStream(&stream, VSL_BRNG_SFMT19937,
		std::chrono::system_clock::now().time_since_epoch().count());
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
			     positions[0].data(), 0, Lx);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
			     positions[1].data(), 0, Ly);
#else
    std::uniform_real_distribution<double> rnd(0, 1.);

	for (long i = 0 ; i < n_parts ; ++i) {
		positions[0][i] = Lx * rnd(rng);
		positions[1][i] = Ly * rnd(rng);
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
	calcForces();

	for (long i = 0 ; i < n_parts ; ++i) {
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

void State::writePos(std::ostream &stream) const {
	for (long i = 0 ; i < n_parts ; ++i) {
		stream << positions[0][i] << " ";
	}
	for (long i = 0 ; i < n_parts ; ++i) {
		stream << positions[1][i] << " ";
	}
	stream << "\n";
}

/*
 * \brief Compute the forces between the particles.
 */
void State::calcForces() {
    for (long i = 0 ; i < n_parts ; ++i) {
		forces[0][i] = 0;
		forces[1][i] = 0;
    }

	calcDists(); // Compute distances
	// Ewald
	//ewald.computeForces(positions[0], positions[1], forces[0], forces[1]);
	ewald.computeForces(positions[0], positions[1], dists[0], dists[1],
			            forces[0], forces[1]);

	// WCA
	/*for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			calcWCAForce(i, j);
		}
	}*/
}

//! Compute internal force between particles i and j (WCA potential)
void State::calcWCAForce(const long i, const long j) {
	double dx = positions[0][i] - positions[0][j];
	double dy = positions[1][i] - positions[1][j];
	// We want the periodized interval to be centered in 0
	pbcSym(dx, Lx);
	pbcSym(dy, Ly);
	double dr2 = (dx * dx + dy * dy) / sigma2;

	if(dr2 * (TWOONETHIRD - dr2) > 0.) {
		double u = WCA_strength / sigma2;
		u *= (48. * pow(dr2, -7.) - 24.*pow(dr2, -4.)); 
		double fx = u * dx;
		double fy = u * dy;

		forces[0][i] += fx;
		forces[0][j] -= fx;
		forces[1][i] += fy;
		forces[1][j] -= fy;
	}
}

// Used only for tests outside of the class
void State::calcDists(const std::vector<double> &pos_x,
		const std::vector<double> &pos_y,
		std::vector<double> &dists_x, std::vector<double> &dists_y,
		double Lx, double Ly) {
	long n_parts = pos_x.size();
	volatile union i_cast u;
	double dx, dy;
	long k = 0;
	dists_x.resize(n_parts * (n_parts - 1) / 2);
	dists_y.resize(n_parts * (n_parts - 1) / 2);

	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			dx = pos_x[i] - pos_x[j];
			dy = pos_y[i] - pos_y[j];
			u.d = (dx / Lx) + 6755399441055744.0;
			dx -= Lx * ((int) u.i[0]);
			u.d = (dy / Ly) + 6755399441055744.0;
			dy -= Ly * ((int) u.i[0]);
			dists_x[k] = dx;
			dists_y[k++] = dy;
		}
	}
}

void State::calcDists() {
	double dx, dy;
	long k = 0;
	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			dx = positions[0][i] - positions[0][j];
			dy = positions[1][i] - positions[1][j];
			enforcePBC(dx, dy);
			dists[0][k] = dx;
			dists[1][k++] = dy;
		}
	}
}

void State::enforcePBC() {
#ifdef USE_MKL
	pbcMKL(positions[0], Lx, aux_x, n_parts);
	pbcMKL(positions[1], Ly, aux_y, n_parts);
#else
	for (long i = 0 ; i < n_parts ; ++i) {
		pbc(positions[0][i], Lx);
		pbc(positions[1][i], Ly);
	}
#endif
}

void State::enforcePBC(double &x, double &y) {
	// dirty but efficient when IEEE 754
	volatile union i_cast u;
	u.d = (fac_x * x) + 6755399441055744.0;
	x -= Lx * ((int) u.i[0]);
	u.d = (fac_y * y) + 6755399441055744.0;
	y -= Ly * ((int) u.i[0]);
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
