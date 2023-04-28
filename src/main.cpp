#include "simul.h"

/*!
 * \brief Main function
 *
 * Create and run the simulation.
 */
int main(int argc, char **argv) {
	Simul simulation(argc, argv);

	if (simulation.getStatus() == SIMUL_INIT_HELP) {
		return 0;
	} else if (simulation.getStatus() == SIMUL_INIT_FAILED) {
		return 1;
	}

	simulation.print(std::cout);
	simulation.run();

	return 0;
}
