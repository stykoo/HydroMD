#include <iostream>
#include <cmath>
#include "ewald.h"

Ewald::Ewald(double _Lx, double _Ly, double _alpha, long _N, bool _verbose) :
	Lx(_Lx), Ly(_Ly), fac_x(1./_Lx), fac_y(1./_Ly), alpha(_alpha), N(_N),
	verbose(_verbose), alpha2(_alpha * _alpha)
{
	double safety = -std::log(EWALD_ERROR);

	// Initializations for real space
	rRange2 = safety / alpha2;
	double rRange = std::sqrt(rRange2);
	hi_x = (long) std::floor(rRange / Lx + 0.5);
	hi_y = (long) std::floor(rRange / Ly + 0.5);
	no_image =  (hi_x == 0 && hi_y == 0);

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

	ones.assign(2 * N, 1.);
#ifdef USE_MKL
	sp = (double *) mkl_malloc(N * nk * sizeof(double), ALIGN);
	Sr = (double *) mkl_malloc(2 * nk * sizeof(double), ALIGN);
	Si = &Sr[nk]; // Both arrays are one after the other
	cc = (double *) mkl_malloc(2 * N * nk * sizeof(double), ALIGN);
	ss = &cc[N * nk]; // Both arrays are one after the other
#else
	Sr = (double *) malloc(nk * sizeof(double));
	Si = (double *) malloc(nk * sizeof(double));
	sp = (double *) malloc(N * nk * sizeof(double));
	cc = (double *) malloc(N * nk * sizeof(double));
	ss = (double *) malloc(N * nk * sizeof(double));
#endif
	

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

Ewald::~Ewald() {
#ifdef USE_MKL
	mkl_free(Sr);
	//mkl_free(Si);
	mkl_free(sp);
	mkl_free(cc);
	//mkl_free(ss);
#else
	free(Sr);
	free(Si);
	free(sp);
	free(cc);
	free(ss);
#endif
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
				force_self_x += (2 * px - 1) / r2 / (2. * M_PI);
				force_self_y += 2 * py / r2 / (2. * M_PI);
			}
		}
	}
}

/* Direct computation of the forces */
void Ewald::computeForcesNaive(
		const std::vector<double> &pos_x, const std::vector<double> &pos_y,
		std::vector<double> &forces_x, std::vector<double> &forces_y,
		double Lx, double Ly) {
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
			if (dx > Lx / 2.)
				dx -= Lx;
			else if (dx < -Lx / 2.)
				dx += Lx;
			if (dy > Ly / 2.)
				dy -= Ly;
			else if (dy < -Ly / 2.)
				dy += Ly;
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
						fx += (2 * px - 1) / r2 / (2. * M_PI);
						fy += 2 * py / r2 / (2. * M_PI);
					}
				}
			}
			forces_x[i] += fx;
			forces_y[i] += fy;
			if (j != i) {
				forces_x[j] += fx;
				forces_y[j] += fy;
			}
		}
	}
}

void Ewald::computeForces(
		const std::vector<double> &pos_x, const std::vector<double> &pos_y,
		const std::vector<double> &dists_x, const std::vector<double> &dists_y,
		std::vector<double> &forces_x, std::vector<double> &forces_y) {
	// Provide distances for speed

	// Self interaction
	for (long i = 0 ; i < N ; ++i) {
		forces_x[i] = force_self_x;
		forces_y[i] = force_self_y;
	}

	// Fourier space contribution
	addFourierForces(pos_x, pos_y, forces_x, forces_y);

	// Real space contribution
	addRealForces(dists_x, dists_y, forces_x, forces_y);
}

void Ewald::addFourierForces(
		const std::vector<double> &pos_x, const std::vector<double> &pos_y,
		std::vector<double> &forces_x, std::vector<double> &forces_y) {

	calcScalarProd(pos_x, pos_y);

	// Sin and cos
#ifdef USE_MKL
	//vdSinCos(N * nk, sp.data(), ss.data(), cc.data());
	vdSinCos(N * nk, sp, ss, cc);
#else
	for (long k = 0 ; k < N * nk ; ++k) {
		cc[k] = std::cos(sp[k]);
		ss[k] = std::sin(sp[k]);
	}
#endif

	calcStructFac();

	// Forces
	double z;
	for (long q = 0 ; q < nk ; ++q) {
		for (long i = 0 ; i < N ; ++i) {
			z = (Sr[q] * cc[N*q+i]) + (Si[q] * ss[N*q+i]) - 1.; // Exclude j=i
			forces_x[i] += z * fCoeffs_x[q];
			forces_y[i] += z * fCoeffs_y[q];
		}
	}
	/*for (long q = 0 ; q < nk ; ++q) {
		for (long i = 0 ; i < N ; ++i) {
			sp[N*q+i] = (Sr[q] * cc[N*q+i]) + (Si[q] * ss[N*q+i]) - 1.;
		}
	}
	cblas_dgemv(CblasRowMajor, CblasTrans, nk, N, 1., sp, N,
		        fCoeffs_x.data(), 1, 1., forces_x.data(), 1);
	cblas_dgemv(CblasRowMajor, CblasTrans, nk, N, 1., sp, N,
		        fCoeffs_y.data(), 1, 1., forces_y.data(), 1);*/

	// Constant contribution
	for (long i = 0 ; i < N ; ++i) {
		forces_x[i] += (N - 1) * fForceAvg_x;
	}
}

void Ewald::calcScalarProd(
		const std::vector<double> &pos_x, const std::vector<double> &pos_y) {
	// Scalar products
	for (long q = 0 ; q < nk ; ++q) {
		for (long i = 0 ; i < N ; ++i) {
			sp[N*q+i] = pos_x[i] * fVecs_x[q] + pos_y[i] * fVecs_y[q];
		}
	}
	/* for (long k = 0 ; k < N * nk ; ++k) {
		sp[k] = 0;
	}
	cblas_dger(CblasRowMajor, nk, N, 1., fVecs_x.data(), 1, pos_x.data(), 1,
			   sp, N);
	cblas_dger(CblasRowMajor, nk, N, 1., fVecs_y.data(), 1, pos_y.data(), 1,
			   sp, N); */
}

void Ewald::calcStructFac() {
	// Structure factor
#if USE_MKL
	// We process both cc/ss and Sr/Si in one go
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*nk, N, 1., cc, N,
		        ones.data(), 1, 0., Sr, 1);
	//cblas_dgemv(CblasRowMajor, CblasNoTrans, nk, N, 1., ss, N,
	//            ones.data(), 1, 0., Si, 1);
#else
	for (long q = 0 ; q < nk ; ++q) {
		Sr[q] = 0.;
		Si[q] = 0.;
	}
	for (long q = 0 ; q < nk ; ++q) {
		for (long i = 0 ; i < N ; ++i) {
			Sr[q] += cc[q*N+i];
			Si[q] += ss[q*N+i];
		}
	}
#endif
}

void Ewald::addRealForces(
		const std::vector<double> &dists_x, const std::vector<double> &dists_y,
		std::vector<double> &forces_x, std::vector<double> &forces_y) {
	double dx, dy, dr2, fx, fy;

	long k = 0;
	for (long i = 0 ; i < N ; ++i) {
		for (long j = 0 ; j < i ; ++j) {
			dx = dists_x[k];
			dy = dists_y[k++];
			dr2 = dx * dx + dy * dy;
			
			if (dr2 < rRange2) {
				if (no_image)
					realForceNoImage(dx, dy, dr2, fx, fy); // without images
				else
					realForce(dx, dy, fx, fy);
				forces_x[i] += fx;
				forces_y[i] += fy;
				forces_x[j] += fx;
				forces_y[j] += fy;
			}
		}
	}
}

void Ewald::realForceNoImage(double dx, double dy, double dr2,
		double &fx, double &fy) {
	// No images, no check for range
	double pref, px, py;
	fx = 0.;
	fy = 0.;

	pref = std::exp(-alpha2 * dr2) / (2. * M_PI);
	px = dx * dx / dr2;
	py = dx * dy / dr2;
	fx += pref * 2 * alpha2 * px;
	fy += pref * 2 * alpha2 * py;
	px = (2 * px - 1) / dr2;
	py = 2 * py / dr2;
	fx += pref * px;
	fy += pref * py;
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
