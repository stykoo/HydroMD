#include <fstream>
#include <exception>
#include <iomanip>
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
		("lx", po::value<double>(&len_x)->default_value(1.0),
		 "Length in x direction")
		("ly", po::value<double>(&len_y)->default_value(1.0),
		 "Length in y direction")
		("parts,n", po::value<long>(&n_parts)->required(),
		 "Number of particles")
		("radius,a", po::value<double>(&radius)->required(),
		 "Particle radius")
		("hydro", po::value<double>(&hydro_strength)->default_value(1.0),
		 "Strength of the hydrodynamic interaction")
		("theta", po::value<double>(&theta_deg)->default_value(0.),
		 "Angle of the hydrodynamic interaction")
		("wca", po::value<double>(&WCA_strength)->default_value(1.0),
		 "Strength of WCA potential")
		("mag", po::value<double>(&mag_strength)->default_value(1.0),
		 "Strength of magnetic interaction")
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
		("extend,e", po::value<std::string>(&extend)->default_value(""),
		 "Name of the file to extend if any")
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
		|| notPositive(WCA_strength, "wca")
		|| notStrPositive(alpha_ew, "alpha")
		|| notStrPositive(dt, "dt")
		|| notPositive(n_iters, "n_iters")) {
		status = SIMUL_INIT_FAILED;
		return;
	}

	std::ifstream infile(output + "_pos.dat");
	if (infile.good()) {
		std::cerr << "Warning: " << output << "_pos.dat already exists. \n"
		          << "Please remove it before running the simulation."
				  << std::endl;
		status = SIMUL_INIT_FAILED;
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
	double theta = theta_deg * M_PI / 180.;
	State state(len_x, len_y, n_parts, radius, hydro_strength, theta,
			    WCA_strength, mag_strength, dt, alpha_ew, extend);

	if (verbose) {
		std::cout << "System initialized\n";
	}
	
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
			  << n_parts << ", radius=" << radius << ", hydro_strength="
			  << hydro_strength << ", theta=" << theta_deg << ", WCA_strength="
			  << WCA_strength << ", mag_strength=" << mag_strength
			  << ", alpha=" << alpha_ew << ", dt=" << dt << ", n_iters="
			  << n_iters << ", n_iters_th=" << n_iters_th << ", skip="
			  << skip << "\n";
	stream << std::endl;
}

int testEwald() {
	double Lx = 1.5, Ly = 1.;
	//double Lx = 1., Ly = 1.;
	double hydro_strength = 1.;
	double theta = 78. * M_PI / 180.;
	std::vector<double> alphas = {1, 1.5, 2., 2.5, 3};

	long N = 5;
	/*std::vector<double> pos_x = {0.2, 0.4};
	std::vector<double> pos_y = {0.3, 0.7};
	std::vector<double> force_x(2, 0.), force_y(2, 0.);*/
	/*std::vector<double> pos_x = {0.77354243, 0.2401132, 0.76776071};
	std::vector<double> pos_y = {0.09785807, 0.9736114, 0.23561104};
	std::vector<double> force_x(3, 0.), force_y(3, 0.);*/
	std::vector<double> pos_x = {0.77354243, 0.2401132, 0.76776071, 0.27041645,
		                         0.74654512};
	std::vector<double> pos_y = {0.09785807, 0.9736114, 0.23561104, 0.83511414,
		                         0.75500363};
	std::vector<double> force_x0(N, 0.), force_y0(N, 0.);
	std::vector<double> force_x(N, 0.), force_y(N, 0.);
	std::vector<double> dists_x, dists_y;

	std::cout << "Testing Ewald method: N=" << N << ", Lx=" << Lx << ", Ly="
		<< Ly << ", hydro_strength=" << hydro_strength << ", theta=" << theta
		<< "\n";

	State::calcDists(pos_x, pos_y, dists_x, dists_y, Lx, Ly);

	Ewald::computeForcesNaive(pos_x, pos_y, force_x0, force_y0, hydro_strength,
			                  theta, Lx, Ly);
	for (size_t i = 0 ; i < pos_x.size() ; ++i) {
		std::cout << std::setprecision(10)
			<< force_x0[i] << " " << force_y0[i] << "\n";
	}

	double max_err = 0., err;
	for (double alpha : alphas) {
		Ewald ew(Lx, Ly, alpha, hydro_strength, theta, N, false);
		ew.computeForces(pos_x, pos_y, dists_x, dists_y, force_x, force_y);
		for (size_t i = 0 ; i < pos_x.size() ; ++i) {
			std::cout << std::setprecision(10)
				<< force_x[i] << " " << force_y[i] << "\n";
			err = std::abs(force_x[i] - force_x0[i]);
			if (err > max_err)
				max_err = err;
			err = std::abs(force_y[i] - force_y0[i]);
			if (err > max_err)
				max_err = err;
		}
	}

	std::cout << "Maximal error: " << max_err << "\n";

	return 0;
}
