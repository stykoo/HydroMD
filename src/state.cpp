#include <cmath>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include "state.h"

/*!
 * \brief Constructor of State
 *
 * Initializes the state of the system: particles randomly placed in a 2d box.
 */
State::State(double _len_x, double _len_y, long _n_parts, double _a,
		     double _hydro_strength, double _theta, double _WCA_strength,
			 double _mag_strength, double _dt, double _alpha_ew,
			 std::string extend) :
	Lx(_len_x), Ly(_len_y), fac_x(1. / _len_x), fac_y(1. / _len_y),
	n_parts(_n_parts), sigma2(4 * _a * _a / TWOONETHIRD),
	WCA_strength(_WCA_strength), mag_strength(_mag_strength), dt(_dt),
	ewald(_len_x, _len_y, _alpha_ew, _hydro_strength, _theta, _n_parts)
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
#endif

	if (extend.empty()) {
		// Generate new configuration
#ifdef USE_MKL
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
					 positions[0].data(), 0, Lx);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_parts,
					 positions[1].data(), 0, Ly);
#else
		std::uniform_real_distribution<double> rnd(0, 1.);

		for (long i = 0 ; i < n_parts ; ++i) {
			positions[0][i] = Lx * rnd(rng);
			positions[1][i] = Ly * rnd(rng);
		}
#endif
		relax(); // Separate particles from one another
	} else {
		// Load positions from file
		loadPos(extend);
	}
}


void State::relax() {
    for (long i = 0 ; i < n_parts ; ++i) {
		forces[0][i] = 0;
		forces[1][i] = 0;
    }
	calcDists();
	double minDist2 = minDistSq();

	long n = 0;
	while (minDist2 < HARMONIC_FAC * TWOONETHIRD * sigma2
		   && n < MAX_RELAX_STEPS) {
		computeHarmonicForces();
		for (long i = 0 ; i < n_parts ; ++i) {
			positions[0][i] += dt * forces[0][i];
			positions[1][i] += dt * forces[1][i];
			forces[0][i] = 0;
			forces[1][i] = 0;
		}
		enforcePBC();
		calcDists();
		minDist2 = minDistSq();
		++n;
	}
	//std::cout << "Min dist: " << sqrt(minDist2) << " (" << n << " iters)\n";
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

void State::loadPos(std::string fname) {
	std::ifstream ifile(fname);
	if (!ifile.is_open()) {
		std::cerr << "Could not find " << fname << "\n";
		return;
	}
	std::string line0, line;
	// Load last line (there should be a faster way)
	while(std::getline(ifile, line0)) {
		line = line0;
	}
	std::istringstream iss(line);

	for (long i = 0 ; i < n_parts ; ++i) {
		iss >> positions[0][i];
	}
	for (long i = 0 ; i < n_parts ; ++i) {
		iss >> positions[1][i];
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

	// Compute distances
	calcDists();

	// Ewald
	ewald.computeForces(positions[0], positions[1], dists[0], dists[1],
			            forces[0], forces[1]);

	// WCA
	computeWCAForces();

	// Magnetic
	computeMagneticForces();
}

void State::computeHarmonicForces() {
	// Assuming distances are already computed
	double dx, dy, dr2, fx, fy, u;
	long k = 0;
	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			dx = dists[0][k];
			dy = dists[1][k++];
			dr2 = (dx * dx + dy * dy) / (TWOONETHIRD * sigma2);

			if(dr2 * (1. - dr2) > 0.) {
				u = HARMONIC_STRENGTH * (1.0 / std::sqrt(dr2) - 1.0);
				fx = u * dx;
				fy = u * dy;

				forces[0][i] += fx;
				forces[0][j] -= fx;
				forces[1][i] += fy;
				forces[1][j] -= fy;
			}
		}
	}
}

void State::computeWCAForces() {
	// Assuming distances are already computed
	double dx, dy, dr2, fx, fy, u;
	long k = 0;
	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			dx = dists[0][k];
			dy = dists[1][k++];
			dr2 = (dx * dx + dy * dy) / sigma2;

			if(dr2 * (TWOONETHIRD - dr2) > 0.) {
				u = WCA_strength / sigma2;
				u *= (48. * pow(dr2, -7.) - 24.*pow(dr2, -4.)); 
				fx = u * dx;
				fy = u * dy;

				forces[0][i] += fx;
				forces[0][j] -= fx;
				forces[1][i] += fy;
				forces[1][j] -= fy;
			}
		}
	}
}

void State::computeMagneticForces() {
	// Assuming distances are already computed
	double dx, dy, dr2, fx, fy, u;
	long k = 0;
	for (long i = 0 ; i < n_parts ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			dx = dists[0][k];
			dy = dists[1][k++];
			dr2 = dx * dx + dy * dy;

			if(dr2 * (Ly * Ly / 4. - dr2) > 0.) {
				u = 3 * mag_strength / pow(dr2, 2.5);
				fx = u * dx;
				fy = u * dy;

				forces[0][i] += fx;
				forces[0][j] -= fx;
				forces[1][i] += fy;
				forces[1][j] -= fy;
			}
		}
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

double State::minDistSq() const {
	double minDist2 = Lx, dr2;
	for (long k = 0 ; k < n_parts * (n_parts - 1) / 2 ; ++k) {
		dr2 = dists[0][k] * dists[0][k] + dists[1][k] * dists[1][k];
		if (dr2 < minDist2)
			minDist2 = dr2;
	}
	return minDist2;
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
