#include <iostream>
#include <cmath>
#include "ewald.h"

Ewald::Ewald(double _Lx, double _Ly, double _alpha, bool _verbose) :
	Lx(_Lx), Ly(_Ly), alpha(_alpha), verbose(_verbose),
	Lx2(_Lx / 2.), Ly2(_Ly / 2), alpha2(_alpha * _alpha)
{
	double safety = -std::log(EWALD_ERROR);

	// Initializations for real space
	rRange2 = safety / alpha2;
	double rRange = std::sqrt(rRange2);
	hi_x = (long) std::ceil(rRange / Lx);
	hi_y = (long) std::ceil(rRange / Ly);

	// Initializations for Fourier space
	double twopi = 2 * M_PI;
	double nf2Max = 4 * alpha2 * safety; // Max norm of Fourier vectors
	long fRangeX = (long) (std::sqrt(nf2Max) * Lx / twopi);
	long fRangeY = (long) (std::sqrt(nf2Max) * Ly / twopi);
	fForceAvg_x = -1. / (2 * Lx * Ly);

	// Building Fourier vectors
	double kx, ky, k2, c;
	for (long mx=-fRangeX ; mx<=fRangeX ; ++mx) {
		kx = 2 * M_PI * mx / Lx;
		for (long my=1 ; my<=fRangeY ; ++my) {
			ky = 2 * M_PI * my / Ly;
			k2 = kx * kx + ky * ky;
			if (k2 < nf2Max) {
				fVecs_x.push_back(kx);
				fVecs_y.push_back(ky);
				c = kx  * std::exp(-k2 / (4. * alpha2)) / k2;
				c *= -2 / (Lx * Ly); // Factor 2 because of symmetry
				fCoeffs_x.push_back(c * kx);
				fCoeffs_y.push_back(c * ky);
			}
		}
	}
	ky = 0.;
	for (long mx=1 ; mx<=fRangeX ; ++mx) {
		kx = 2 * M_PI * mx / Lx;
		k2 = kx * kx;
		if (k2 < nf2Max) {
			fVecs_x.push_back(kx);
			fVecs_y.push_back(ky);
			c = kx  * std::exp(-k2 / (4. * alpha2)) / k2;
			c *= -2 / (Lx * Ly);
			fCoeffs_x.push_back(c * kx);
			fCoeffs_y.push_back(c * ky);
		}
	}
	nk = fVecs_x.size();
	// Structure factor (real and imaginary parts)
	Sr.resize(nk);
	Si.resize(nk);
	

	// Self interaction
	computeSelfInteraction();

	std::cout << "Ewald: Lx=" << Lx << ", Ly=" << Ly << ", alpha=" << alpha
		<< ", rRange=" << rRange << ", fRangeX=" << fRangeX << ", fRangeY="
		<< fRangeY << "\n";

	if (verbose) {
		std::cout << "force_self = " << force_self_x << ", " << force_self_y
			<< "\n\n";
		std::cout << nk << " Fourier vectors" << "\n";
		for (long i = 0 ; i < nk ; ++i) {
			std::cout << fVecs_x[i] << ", " << fVecs_y[i] << "\n";
		}
		std::cout << "\n";
		std::cout << nk << " Fourier coefficients" << "\n";
		for (long i = 0 ; i < nk ; ++i) {
			std::cout << fCoeffs_x[i] << ", " << fCoeffs_y[i] << "\n";
		}
		std::cout << "\n";
	}
}

void Ewald::computeSelfInteraction() {
	// Direct computation of the interaction of a particle with itself
	double rRange2A = 1. / EWALD_ERROR;
	double rRangeA = sqrt(rRange2A);
	long hi_xA = (long) std::ceil(rRangeA / Lx);
	long hi_yA = (long) std::ceil(rRangeA / Ly);
	//std::cout << rRangeA << ", " << hi_xA << ", " << hi_yA << "\n";
	double x, y, r2, px, py; 

	force_self_x = 0.;
	force_self_y = 0.;

	for (long a = -hi_xA ; a <= hi_xA ; ++a) {
		x = a * Lx;
		for (long b = -hi_yA ; b <= hi_yA ; ++b) {
			y = b * Ly;
			r2 = x * x + y * y;
			if (r2 > 0 && r2 < rRange2A) {
				px = x * x / r2;
				py = x * y / r2;
				force_self_x += (2 * px - 1) / r2;
				force_self_y += 2 * py / r2;
			}
		}
	}
}

/* Direct computation of the forces */
void Ewald::computeForcesNaive(
		const std::vector<double> &pos_x, const std::vector<double> &pos_y,
		std::vector<double> &forces_x, std::vector<double> &forces_y) {
	long N = pos_x.size();

	for (long i = 0 ; i < N ; ++i) {
		forces_x[i] = 0;
		forces_y[i] = 0;
	}

	double rRange2A = 1. / EWALD_ERROR;
	double rRangeA = sqrt(rRange2A);
	long hi_xA = (long) std::ceil(rRangeA / Lx);
	long hi_yA = (long) std::ceil(rRangeA / Ly);
	double dx, dy, x, y, r2, px, py, fx, fy; 

	for (long i = 0 ; i < N ; ++i) {
		for (long j = 0 ; j <= i ; ++j) {
			dx = pos_x[i] - pos_x[j];
			dy = pos_y[i] - pos_y[j];
			enforcePBC(dx, dy);
			fx = 0;
			fy = 0;
			for (long a = -hi_xA ; a <= hi_xA ; ++a) {
				x = dx - a * Lx;
				for (long b = -hi_yA ; b <= hi_yA ; ++b) {
					y = dy - b * Ly;
					r2 = x * x + y * y;
					if (r2 > 0 && r2 < rRange2A) {
						px = x * x / r2;
						py = x * y / r2;
						fx += (2 * px - 1) / r2;
						fy += 2 * py / r2;
					}
				}
			}
			forces_x[i] += fx;
			forces_y[i] += fy;
			forces_x[j] += fx;
			forces_y[j] += fy;
		}
	}
}

void Ewald::computeForces(
		const std::vector<double> &pos_x, const std::vector<double> &pos_y,
		std::vector<double> &forces_x, std::vector<double> &forces_y) {
	long N = pos_x.size();

	// Self interaction
	for (long i = 0 ; i < N ; ++i) {
		forces_x[i] = force_self_x;
		forces_y[i] = force_self_y;
	}

	// Fourier space contribution
	addFourierForces(pos_x, pos_y, forces_x, forces_y);

	// Real space contribution
	addRealForces(pos_x, pos_y, forces_x, forces_y);
}


void Ewald::addFourierForces(
		const std::vector<double> &pos_x, const std::vector<double> &pos_y,
		std::vector<double> &forces_x, std::vector<double> &forces_y) {
	long N = pos_x.size();
	long Nall = N * nk; 
	sp.resize(Nall); // scalar products G.R
	cc.resize(Nall); // cos(G.r)
	ss.resize(Nall); // sin(G.r)

	// Scalar products
	long k = 0;
	for (long q = 0 ; q < nk ; ++q) {
		for (long i = 0 ; i < N ; ++i) {
			sp[k] = pos_x[i] * fVecs_x[q] + pos_y[i] * fVecs_y[q];
			++k;
		}
	}

	// Sin and cos
#ifdef USE_MKL
	vdSinCos(Nall, sp.data(), ss.data(), cc.data());
#else
	for (k = 0 ; k < Nall ; ++k) {
		cc[k] = std::cos(sp[k]);
		ss[k] = std::sin(sp[k]);
	}
#endif

	// Structure factor
	k = 0;
	for (long q = 0 ; q < nk ; ++q) {
		Sr[q] = 0.;
		Si[q] = 0.;
		for (long i = 0 ; i < N ; ++i) {
			Sr[q] += cc[k];
			Si[q] += ss[k];
			++k;
		}
	}
	
	// Forces
	double z;
	k = 0;
	for (long q = 0 ; q < nk ; ++q) {
		for (long i = 0 ; i < N ; ++i) {
			z = Sr[q] * cc[k] + Si[q] * ss[k];
			z -= 1.; // Exclude j = i term
			forces_x[i] += z * fCoeffs_x[q];
			forces_y[i] += z * fCoeffs_y[q];
			++k;
		}
	}

	// Constant contribution
	for (long i = 0 ; i < N ; ++i) {
		forces_x[i] += fForceAvg_x;
	}
}

void Ewald::addRealForces(
		const std::vector<double> &pos_x, const std::vector<double> &pos_y,
		std::vector<double> &forces_x, std::vector<double> &forces_y) {
	long N = pos_x.size();
	double dx, dy, fx, fy;

	for (long i = 0 ; i < N ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			dx = pos_x[i] - pos_x[j];
			dy = pos_y[i] - pos_y[j];
			enforcePBC(dx, dy);
			realForce(dx, dy, fx, fy);
			forces_x[i] += fx;
			forces_y[i] += fy;
			forces_x[j] += fx;
			forces_y[j] += fy;
		}
	}
}

void Ewald::enforcePBC(double &x, double &y) {
	if (x > Lx2)
		x -= Lx;
	else if (x < -Lx2)
		x += Lx;
	if (y > Ly2)
		y -= Ly;
	else if (y < -Ly2)
		y += Ly;
}

void Ewald::realForce(double dx, double dy, double &fx, double &fy) {
	double dx0, dy0, dr2, pref, px, py;
	fx = 0.;
	fy = 0.;

	// Sum over periodic boxes
	for (long a = -hi_x ; a <= hi_x ; ++a) {
		dx0 = dx - a * Lx;
		for (long b = -hi_y ; b <= hi_y ; ++b) {
			dy0 = dy - b * Ly;
			dr2 = dx0 * dx0 + dy0 * dy0;
			if (dr2 > 0 && dr2 < rRange2) {
				pref = std::exp(-alpha2 * dr2) / (2. * M_PI);
				px = dx0 * dx0 / dr2;
				py = dx0 * dy0 / dr2;
				fx += pref * 2 * alpha2 * px;
				fy += pref * 2 * alpha2 * py;
				px = (2 * px - 1) / dr2;
				py = 2 * py / dr2;
				fx += pref * px;
				fy += pref * py;
			}
		}
	}
}

int testEwald() {
	double Lx = 1., Ly = 1.;
	double alpha = 3;
	Ewald ew(Lx, Ly, alpha, false);

	std::vector<double> pos_x = {0.2, 0.4};
	std::vector<double> pos_y = {0.3, 0.7};
	std::vector<double> force_x(2, 0.), force_y(2, 0.);

	ew.computeForcesNaive(pos_x, pos_y, force_x, force_y);
	for (long i = 0 ; i < 2 ; ++i) {
		std::cout << force_x[i] << " " << force_y[i] << "\n";
	}
	std::cout << "\n";

	ew.computeForces(pos_x, pos_y, force_x, force_y);
	for (long i = 0 ; i < 2 ; ++i) {
		std::cout << force_x[i] << " " << force_y[i] << "\n";
	}

	return 0;
}
