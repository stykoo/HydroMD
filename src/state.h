#ifndef HYDROMD_STATE_H_
#define HYDROMD_STATE_H_

#include <vector>
#include <array>
#include <iostream>
#include "ewald.h"

#ifdef USE_MKL
	#include "mkl.h"
	#include "mkl_vsl.h"
#else
	#include <random>
#endif

// 2^(1/6)
// #define TWOONESIXTH 1.12246204830937298143 
// 2(1/3)
#define TWOONETHIRD 1.25992104989487316477
#define HARMONIC_STRENGTH 10
#define HARMONIC_FAC 0.9
#define MAX_RELAX_STEPS 20000


/*!
 * \brief Class for the state of the system
 *
 * This class takes care of the initialization
 * and the evolution of the state of the system.
 */
class State {
	public:
		//! Constructor of State
		State(double _len_x, double _len_y, long _n_parts, double _a,
		      double _hydro_strength, double _theta, double _WCA_strength,
		      double _mag_strength, double _dt, double _alpha_ew,
		      std::string extend);
		~State() {
#ifdef USE_MKL
			vslDeleteStream(&stream);
#endif
		}
		void evolve(); //!< Do one time step

		//! Get the x coordinate of the positions 
		const std::vector<double> & getPosX() const {
			return positions[0];
		}
		//! Get the y coordinate of the positions 
		const std::vector<double> & getPosY() const {
			return positions[1];
		}

		void dump() const; //!< Dump the positions
		void dumpHydroForces(); //!< Dump the hydrodynamic forces
		void writePos(std::ostream &stream) const;
		//! For use outside of the class
		static void calcDists(const std::vector<double> &pos_x,
				const std::vector<double> &pos_y,
				std::vector<double> &dists_x, std::vector<double> &dists_y,
				double Lx, double Ly);


	private:
		void relax(); //!< Relax the initial system
		void calcForces(); //!< Compute internal forces
		void calcDists(); //!< Compute distances
		void computeHarmonicForces(); //!< Harmonic repulsive forces
		void computeWCAForces(); //!< WCA repulsive forces
		void computeMagneticForces(); //!< Magnetic forces
		void loadPos(std::string fname);
		void enforcePBC(); //!< Enforce periodic boundary conditions
		void enforcePBC(double &x, double &y); //!< Same on specific numbers
		double minDistSq() const;

		const double Lx, Ly; //!< Length and width of the box
		const double fac_x, fac_y; //!< Inverse length and width of the box
		const long n_parts; //!< Number of particles
		const double sigma2; //!< Square of particle diameter
		const double WCA_strength; //!< Strength of the WCA potential
		const double mag_strength; //!< Strength of the magnetic interaction
		const double dt; //!< Timestep
		Ewald ewald;

#ifdef USE_MKL
		VSLStreamStatePtr stream;
		std::vector<double> aux_x, aux_y;
#else
		std::mt19937 rng; //!< Random number generator
#endif

		//! Positions of the particles
		std::array<std::vector<double>, 2> positions;
		std::array<std::vector<double>, 2> forces;  //!< Internal forces
		std::array<std::vector<double>, 2> dists;
};

// Trick to avoid round ASSUMING LITTLE ENDIAN
union i_cast {double d; int i[2];};
#define double2int(i, d, t)  \
    {volatile union i_cast u; u.d = (d) + 6755399441055744.0; \
    (i) = (t)u.i[0];}


template<typename T>
void pbc(T &x, const T L) {
	x -= L * std::floor(x / L);
}

template<typename T>
inline void pbcSym(T &x, const T L) {
	x -= L * std::round(x / L);
}

#ifdef USE_MKL
void pbcMKL(std::vector<double> &v, const double L, std::vector<double> &aux,
	        const long N);
void pbcSymMKL(std::vector<double> &v, const double L,
		       std::vector<double> &aux, const long N);
#endif


#endif // HYDROMD_STATE_H_
