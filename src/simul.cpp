#include <fstream>
#include <exception>
#include <boost/program_options.hpp>
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
		("radius,a", po::value<double>(&radius)->required(),
		 "Particle radius")
		("eps", po::value<double>(&WCA_strength)->default_value(1.0),
		 "Strength of WCA potential")
		("alpha", po::value<double>(&alpha_ew)->required(),
		 "Parameter of Ewald algorithm")
		("dt,t", po::value<double>(&dt)->required(), "Timestep")
		("iters,I", po::value<long>(&n_iters)->required(),
		 "Number of time iterations")
		("itersTh,J", po::value<long>(&n_iters_th)->default_value(0),
		 "Number of time iterations of thermalization")
		("skip,S", po::value<long>(&skip)->default_value(100),
		 "Iterations between two computations of observables")
		("output,O",
		 po::value<std::string>(&output)->default_value("observables"),
		 "Name of the output file")
		("pos", po::bool_switch(&export_pos), "Export positions")
		("test", po::bool_switch(&test), "Test mode")
		("verbose,v", po::bool_switch(&verbose), "Verbose mode")
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
		|| notStrPositive(radius, "a")
		|| notPositive(WCA_strength, "eps")
		|| notStrPositive(alpha_ew, "alpha")
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
	// Test of Ewald method
	if (test) {
		testEwald();
		return;
	}

	// Initialize the state of the system
	State state(len_x, len_y, n_parts, radius, WCA_strength, dt, alpha_ew);
	
	// Thermalization
	for (long t = 0 ; t < n_iters_th ; ++t) {
		if (verbose) 
			std::cout << t << "\r";
		state.evolve();
	}

	std::ofstream ofile;
	if (export_pos) {
		ofile.open(output + "_pos.dat");
		print(ofile);
		state.writePos(ofile);
	}

	// Time evolution
	for (long t = 0 ; t < n_iters ; ++t) {
		if (verbose) 
			std::cout << t << "\r";
		state.evolve();
		if (t % skip == 0) {
			if (export_pos) {
				state.writePos(ofile);
			}
		}
	}
	if (export_pos) {
		ofile.close();
	}
}

/*!
 * \brief Print the parameters of the simulation
 */
void Simul::print(std::ostream &stream) const {
	stream << "# ";
	stream << "len_x=" << len_x << ", len_y=" << len_y << ", n_parts="
			  << n_parts << ", radius=" << radius <<  ", WCA_strength="
			  << WCA_strength << ", alpha=" << alpha_ew << ", dt=" << dt
			  << ", n_iters=" << n_iters << ", n_iters_th=" << n_iters_th
			  << ", skip=" << skip << "\n";
	stream << std::endl;
}
