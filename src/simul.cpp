#include <exception>
#include <boost/program_options.hpp>
// #include "observables.h"
#include "simul.h"
#include "ewald.h"
#include "state.h"

namespace po = boost::program_options;

/*!
 * \brief Constructor of Simul
 *
 * Initializes the parameters of structure Simul
 * from the command-line arguments using boost::program_options.
 *
 * \param argc Number of arguments
 * \param argv Arguments
 */
Simul::Simul(int argc, char **argv) {
	status = SIMUL_INIT_SUCCESS;

	po::options_description opts("Options");
	opts.add_options()
		("lx", po::value<double>(&len_x)->default_value(1.0), "")
		("ly", po::value<double>(&len_y)->default_value(1.0), "")
		("parts,n", po::value<long>(&n_parts)->required(),
		 "Number of particles")
		("eps,e", po::value<double>(&pot_strength)->default_value(1.0),
		 "Strength of interparticle potential")
		("dt,t", po::value<double>(&dt)->required(), "Timestep")
		("iters,I", po::value<long>(&n_iters)->required(),
		 "Number of time iterations")
		("itersTh,J", po::value<long>(&n_iters_th)->default_value(0),
		 "Number of time iterations of thermalization")
		("skip,S", po::value<long>(&skip)->default_value(100),
		 "Iterations between two computations of observables")
		("output,O",
		 po::value<std::string>(&output)->default_value("observables.h5"),
		 "Name of the output file")
		("test", po::bool_switch(&test), "Test mode")
		("help,h", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts), vars);

		// Display help and exit
		if (vars.count("help")) {
			std::cout << "Usage: " << argv[0] << " options\n";
			std::cout << opts << std::endl;
			status = SIMUL_INIT_HELP;
			return;
		}

        po::notify(vars);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		status = SIMUL_INIT_FAILED;
		return;
	}

	// Check if the values of the parameters are allowed
	if (notStrPositive(len_x, "len_x") || notStrPositive(len_y, "len_y")
		|| notStrPositive(n_parts, "n_parts")
		|| notPositive(pot_strength, "eps")
		|| notStrPositive(dt, "dt")
		|| notPositive(n_iters, "n_iters")) {
		status = SIMUL_INIT_FAILED;
		return;
	}
}

/*!
 * \brief Run the simulation
 *
 * Construct the state of the system and update it for the number
 * of iterations wanted. Also take care of launching the thread for
 * visualization.
 */
void Simul::run() {
	if (status != SIMUL_INIT_SUCCESS) {
		std::cerr << "You should not be runing a failed simulation..."
		          << std::endl;
		return;
	}
	if (test) {
		testEwald();
		return;
	}

	// Initialize the state of the system
	State state(len_x, len_y, n_parts, pot_strength, dt);
	/*Observables obs(len, n_parts, step_r, n_div_angle, less_obs,
					cartesian);*/
	
	// Thermalization
	for (long t = 0 ; t < n_iters_th ; ++t) {
		state.evolve();
	}
	// Time evolution
	for (long t = 0 ; t < n_iters ; ++t) {
		state.evolve();
		/*if (t % skip == 0) {
			obs.compute(&state);
		}*/
	}

	/*obs.writeH5(output, rho, n_parts, pot_strength, temperature, rot_dif,
				activity, dt, n_iters, n_iters_th, skip);*/
	//state.dump();
}

/*!
 * \brief Print the parameters of the simulation
 */
void Simul::print() const {
	std::cout << "# ";
	std::cout << "len_x=" << len_x << ", len_y=" << len_y << ", n_parts="
			  << n_parts << ", pot_strength=" << pot_strength
			  << ", dt=" << dt << ", n_iters=" << n_iters
			  << ", n_iters_th=" << n_iters_th << ", skip=" << skip << "\n";
	std::cout << std::endl;
}
