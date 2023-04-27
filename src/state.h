#ifndef HYDROMD_STATE_H_
#define HYDROMD_STATE_H_

#include <vector>
#include <array>
#include "ewald.h"

#ifdef USE_MKL
	#include "mkl.h"
	#include "mkl_vsl.h"
#else
	#include <random>
#endif

// 2^(1/6)
#define TWOONESIXTH 1.12246204830937298143 


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
		      double _pot_strength, double _dt, double _alpha_ew);
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


	private:
		void calcForces(); //!< Compute internal forces
		void calcWCAForce(const long i, const long j);
		void enforcePBC(); //!< Enforce periodic boundary conditions

		const double len_x, len_y; //!< Length and width of the box
		const long n_parts; //!< Number of particles
		const double a; //!< Particle radius
		const double pot_strength; //!< Strength of the interparticle potential
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
};

/*! 
 * \brief Periodic boundary conditions on a segment
 * 
 * Update x to be between 0 and L.
 *
 * \param x Value
 * \param L Length of the box
 */
template<typename T>
void pbc(T &x, const T L) {
	x -= L * std::floor(x / L);
}

/*! 
 * \brief Periodic boundary conditions on a segment (symmetric)
 * 
 * Update x to be between -L/2 and L/2.
 *
 * \param x Value
 * \param L Length of the box
 */
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
