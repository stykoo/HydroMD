#ifndef HYDROMD_SIMUL_H_
#define HYDROMD_SIMUL_H_

#include <string>
#include <iostream>

//! State of the simulation after initialization
enum SimulInitStatus {
	SIMUL_INIT_SUCCESS, //!< Successful initialization
	SIMUL_INIT_HELP, //!< Display help
	SIMUL_INIT_FAILED //!< Failed initialization
};

/*!
 * \brief Class for simulation
 *
 * This class takes care of both the initialization
 * and the time evolution of the system.
 */
class Simul {
	public:
		Simul(int argc, char **argv); //!< Constructor from arguments
		void run(); //!< Run the simulation
		void print(std::ostream &stream) const; //!< Print the parameters

		//! Get initialization status
		SimulInitStatus getStatus() const { return status; }

	private:
		double len_x; //!< Length in x direction
		double len_y; //!< Length in y direction
		long n_parts; //!< Number of particles
		double radius; //!< Radius of particle
		double WCA_strength; //!< Strength of the WCA potential
		double alpha_ew; //!< Parameter of Ewald algorithm
		double dt; //!< Timestep
		long n_iters; //!< Number of time iterations
		long n_iters_th; //!< Number of time iterations of thermalization
		long skip; //!< Iterations between two computation of observables
		std::string output; //!< Name of the output file
		bool export_pos; //!< Export positions
		bool test; //!< Test mode
		bool verbose; //!< Verbose mode

		SimulInitStatus status; //!< Status after initialization
};

/*!
 * \brief Check if variable is positive or null.
 *
 * Returns true and prints an error message if the variable
 * is not positive.
 *
 * \param a Variable
 * \param name Name of variable
 * \return false if the variable is positive, true otherwise
 */
template<typename T>
bool notPositive(const T &a, std::string name) {
	if (a < T(0)) {
		std::cerr << "Error: " << name << " should be positive."
		          << std::endl;
		return true;
	}
	return false;
}

/*!
 * \brief Check if variable is strictly positive.
 *
 * Returns true and prints an error message if the variable
 * is not positive or is null.
 *
 * \param a Variable
 * \param name Name of variable
 * \return false if the variable is strictly positive, true otherwise
 */
template<typename T>
bool notStrPositive(const T &a, std::string name) {
	if (a <= T(0)) {
		std::cerr << "Error: " << name << " should be strictly positive."
		          << std::endl;
		return true;
	}
	return false;
}

int testEwald();

#endif // HYDROMD_SIMUL_H_
