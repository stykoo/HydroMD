#ifndef HYDROMD_EWALD_H_
#define HYDROMD_EWALD_H_

#include <vector>
#ifdef USE_MKL
	#include "mkl.h"
	#include "mkl_vsl.h"
#endif

#define EWALD_ERROR 1e-7
#define ALIGN 256

class Ewald {
	public:
		Ewald(double _Lx, double _Ly, double _alpha, double _hydro_strength,
			  double _theta, long _N, bool _verbose=false);
		~Ewald();

		static void computeForcesNaive(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y,
			std::vector<double> &forces_x, std::vector<double> &forces_y,
			double hydro_strength, double theta, double Lx, double Ly); 
		void computeForces(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y,
			const std::vector<double> &dists_x,
			const std::vector<double> &dists_y,
			std::vector<double> &forces_x, std::vector<double> &forces_y);

	private:
		void computeSelfInteraction();
		void addRealForces(
			const std::vector<double> &dists_x,
			const std::vector<double> &dists_y,
			std::vector<double> &forces_x, std::vector<double> &forces_y);
		void addFourierForces(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y,
			std::vector<double> &forces_x, std::vector<double> &forces_y);

		void calcScalarProd(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y);
		void calcStructFac();
		void calcDists(
			const std::vector<double> &pos_x, const std::vector<double> &pos_y);
		void realForceNoImage(double dx, double dy, double dr2,
				              double &fx, double &fy);
		void realForce(double dx, double dy, double &fx, double &fy);

		const double Lx, Ly; //!< Dimensions of the box
		const double fac_x, fac_y; //!< Inverse dimensions of the box
		const double alpha; //!< Parameter of Ewald algorithm
		const double hydro_strength; //!< Strength of hydrodynamic interaction
		const double ct, st; //!< cos and sin of the dipole angle
		const long N; //!< Number of particles
		const bool verbose; //!< Verbose mode
		const double alpha2;


		double rRange2; //!< Maximal square distance between particles
		long hi_x, hi_y; //!< Maximal images in x/y direction
		bool no_image; //!< Don't need images

		double fForceAvg_x, fForceAvg_y; //!< Constant contribution
		std::vector<double> fVecs_x, fVecs_y; //!< Reciprocal space vectors G
		//< Fourier space coefficients: (-1/V) (Gx/G^2) exp(-G^2/(4 alpha^2)) G
		std::vector<double> fCoeffs_x, fCoeffs_y;
		long nk; //!< Number of reciprocal space vectors

		double force_self_x, force_self_y; // Force of a particle on itself

		std::vector<double> ones;
		double *sp, *Sr, *Si, *cc, *ss; // C-style array for MKL operations
};

#endif // HYDROMD_EWALD_H_
