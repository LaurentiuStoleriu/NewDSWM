#pragma once
#define npart 1
#define nvar (2*npart)
#define cond_sw /*torq_mod>eps */  t < t_stop

const double Pi = 3.1415926535897932384626433832795;
const double Pis2 = 1.5707963267948966192313216916398;
const double Pix2 = 6.283185307179586476925286766559;
const double miu0 = 1.256637062121919191919191919*1.0e-6;
const double gammap = miu0 * 1.760859708*1.0e+11;

const double alpha = 0.01;
const double Ms_mod = 1.0 * 795774.7154594767;
const double K1_mod = 1.0 * 1.0e+5;
const double R_mod = 1.0e-9;
const double V_mod = 4.0*Pi / 3.0 * pow(8.0*R_mod, 3.0);

const double J_dip = 0.0;
const int n_max_vec = npart;

const double Hk_mod = 2.0*K1_mod / miu0 / Ms_mod;

const double step_tSW = 1.001;//ps
const double eps = 1.0e-5;

const double J_exch = 0.0;
const int n_max_vec_exch = 4;
const double prag_vecini_exch = 0.0;// / 4.0 / Pi / R_mod / R_mod / R_mod;

class Camp {
public:
	double Hamp;
	double theta_h, phi_h;
	double Hx, Hy, Hz;
};
class Particula {
public:
	int index, strat;
	double X, Y, Z, R, L_j;
	double K1, K2, K3;
	double Mx, My, Mz;
	double m, V, Ms;
	double KV;
	double Hk, Hk2, T_Larmor_Hk;
	double T_relax;
	double theta_m, phi_m;
	double theta_ea, phi_ea;
	double eax, eay, eaz;
	double utheax, utheay, utheaz;
	double upheax, upheay, upheaz;
};
class Interactiuni {
public:
	int vecin;
	double coef, xx, xy, xz, yx, yy, yz, zx, zy, zz;
};
class Moment {
public:
	double Mamp;
	double theta_M, phi_M;
	double Mx, My, Mz;
};
struct SW_function_params {
public:
	double Hk_fcn, H_fcn, theta0_fcn;
};

double gtr(double);
double rtg(double);
double SafeAcos(double);
double timp_max(double, double, double);
void camp(double, double, double, Camp*);
void camp_sir(int, double, double, double, std::vector<double>*);
void position_coeficients(Particula, Particula, Interactiuni*);
double uniform_rand(double low, double high);
double normal_rand(double mean, double sigma);
