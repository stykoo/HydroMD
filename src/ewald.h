#ifndef HYDROMD_EWALD_H_
#define HYDROMD_EWALD_H_

#include <vector>
#ifdef USE_MKL
	#include "mkl.h"
	#include "mkl_vsl.h"
#endif

#define EWALD_ERROR 1e-7

class Ewald {
	public:
		Ewald(double _Lx, double _Ly, double _alpha, bool _verbose=false);
		void computeForcesNaive(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y,
			std::vector<double> &forces_x, std::vector<double> &forces_y);
		void computeForces(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y,
			std::vector<double> &forces_x, std::vector<double> &forces_y);

	private:
		void computeSelfInteraction();
		void addRealForces(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y,
			std::vector<double> &forces_x, std::vector<double> &forces_y);
		void addFourierForces(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y,
			std::vector<double> &forces_x, std::vector<double> &forces_y);
		void enforcePBC(double &x, double &y);
		void realForce(double dx, double dy, double &fx, double &fy);

		const double Lx, Ly; //!< Dimensions of the box
		const double alpha; //!< Parameter of Ewald algorithm
		const bool verbose; //!< Verbose mode
		const double Lx2, Ly2, alpha2;


		double rRange2; //!< Maximal square distance between particles
		long hi_x, hi_y; //!< Maximal images in x/y direction

		double fForceAvg_x; //!< Constant contribution
		std::vector<double> fVecs_x, fVecs_y; //!< Reciprocal space vectors G
		//< Fourier space coefficients: (-1/V) (Gx/G^2) exp(-G^2/(4 alpha^2)) G
		std::vector<double> fCoeffs_x, fCoeffs_y;
		long nk; //!< Number of reciprocal space vectors

		double force_self_x, force_self_y; // Force of a particle on itself

		std::vector<double> Sr, Si;
		std::vector<double> sp, cc, ss;
};

int testEwald();

#endif // HYDROMD_EWALD_H_
