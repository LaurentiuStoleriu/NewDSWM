#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "Header_DSWM_MHL.h"

//using namespace std;

double gtr(double x) {
	return x * Pi / 180.0;
}
double rtg(double x) {
	return x * 180.0 / Pi;
}
double uniform_rand(double low, double high) {
	double range = high - low;
	return low + (double)(range * rand() / RAND_MAX);
}
double normal_rand(double mean, double sigma) {
	double u1 = uniform_rand(0.0, 1.0);
	double u2 = uniform_rand(0.0, 1.0);
	double z0 = sqrt(-2.0*log(u1))*cos(Pix2*u2);
	return z0 * sigma + mean;
}

double SafeAcos(double x) {
	if (x < -1.0) x = -1.0;
	else if (x > 1.0) x = 1.0;
	return acos(x);
}
double timp_max(double A, double theta_in, double theta_tar) {

	double p1, p2, p3, p4, Y0, ampl, expo, psi, zeta;
	double th_in = 0, theta_in_c, tmc, u, x, timp_polynom = 0;

	p1 = 1.56 + 1.63*pow(0.15, A);
	p2 = (1.17 + 1.07 / A * exp(-log(A / 2.7)*log(A / 2.7)))*theta_tar;
	p3 = (-0.15 - 6.02 / (2 * A*A - 5 * A + 13.88))*theta_tar*theta_tar;
	p4 = (0.03 + 1.25 / (2 * A*A - 5 * A + 13.7))*theta_tar*theta_tar*theta_tar;
	theta_in_c = p1 + p2 + p3 + p4;

	Y0 = 18.9 + 126.0*pow(0.62, A);
	ampl = 7.27 + 831*pow(0.22, A);
	expo = exp(-0.5*(theta_tar - Pis2)*(theta_tar - Pis2) / 0.66 / 0.66);
	tmc = Y0 + ampl * expo;

	psi = 0.86-0.23*pow(0.59, A);
	zeta = 0.04 + 0.06*A;

	u = 2 * theta_in_c - theta_tar;
	if (theta_in < theta_tar)
		th_in = 2 * theta_tar - theta_in;
	if (theta_in > u)
		th_in = 2 * u - theta_in;
	if (theta_in > theta_tar && theta_in < u)
		th_in = theta_in;

	x = th_in - theta_in_c;
	timp_polynom = tan(psi*x) / tan(zeta);

	return timp_polynom + tmc;
}

void camp(double th, double ph, double a, Camp* F) {
	if (a >= 0.0) {
		F->Hamp = a;
		F->theta_h = gtr(th);
		F->phi_h = gtr(ph);
	}
	if (a < 0.0) {
		F->Hamp = -1.0*a;
		F->theta_h = gtr(180.0 - th);
		if (ph < 180.0)
			F->phi_h = gtr(180.0 + ph);
		if (ph > 180.0)
			F->phi_h = gtr(ph - 180.0);
		if (ph ==180.0)
			F->phi_h = gtr(180.0);
	}

	F->Hx = F->Hamp*sin(F->theta_h)*cos(F->phi_h);
	F->Hy = F->Hamp*sin(F->theta_h)*sin(F->phi_h);
	F->Hz = F->Hamp*cos(F->theta_h);
}
void camp_sir(int nr_frc, double a, double t, double T, std::vector<double> *sir) {
	
	
	/*Aceasta functie genereaza un sir de campuri in A/m
	si le stocheaza intr-un vector*/
	
	
	std::vector<double> F;

	double step_tp = t * 1.0e-12;
	double period = T * 1.0e-12;
	int i_max = (int) (T / t);
	double f = 1.0 / period;



	double h;
	for (int i = 0; i < i_max/2; i++) {
		h = a * (1.0 - 2.0*(double)i / ((double)i_max/2 - 1));
		//h = a * cos(Pix2*i*step_tp*f);
		F.push_back(h);
	}
	for (int i = i_max/2; i >= 0; i--) {
		h = a * (1.0 - 2.0*(double)i / ((double)i_max/2- 1));
		//h = a * cos(Pix2*i*step_tp*f);
		F.push_back(h);
	}
	
	int i_stop;
	for (int index_forc = 1; index_forc <= nr_frc; index_forc++) {

		i_stop = i_max/4+(index_forc)*25;
		for (int i=1; i < i_stop; i++) {
			h = a * (1.0 - 2.0*(double)i / ((double)i_max/2 - 1));
			F.push_back(h);
		}
		for (int i = i_stop; i >= 0; i--) {
			h = a * (1.0 - 2.0*(double)i / ((double)i_max/2 - 1));
			F.push_back(h);
		}
	}
	*sir = F;
}

void position_coeficients(Particula P1, Particula P2, Interactiuni *Pos_Coef) {
	double dist = sqrt((P2.X - P1.X)*(P2.X - P1.X) + (P2.Y - P1.Y)*(P2.Y - P1.Y) + (P2.Z - P1.Z)*(P2.Z - P1.Z));
	double rx = (P2.X - P1.X) / dist;
	double ry = (P2.Y - P1.Y) / dist;
	double rz = (P2.Z - P1.Z) / dist;
	Pos_Coef->coef = 1.0 / 4.0 / Pi / dist / dist / dist;
	Pos_Coef->xx = 3.0 * rx * rx - 1.0;
	Pos_Coef->xy = 3.0 * rx * ry;
	Pos_Coef->xz = 3.0 * rx * rz;
	Pos_Coef->yx = Pos_Coef->xy;
	Pos_Coef->yy = 3.0 * ry * ry - 1.0;
	Pos_Coef->yz = 3.0 * ry * rz;
	Pos_Coef->zx = Pos_Coef->xz;
	Pos_Coef->zy = Pos_Coef->yz;
	Pos_Coef->zz = 3.0 * rz * rz - 1.0;
}
