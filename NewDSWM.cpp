#include <iostream>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#include <string>

/*#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>*/

#include "Header_DSWM_MHL.h"

using std::cout; using std::cin; using std::endl;
using std::ofstream;
using std::vector; using std::complex;

typedef vector<double> doublev;
typedef complex<double> doublec;

double torq_mod = 1.0;

Camp H;
doublev sir_camp;

Particula P[npart], P_init[npart];

Interactiuni Pozitie[npart][n_max_vec];
Interactiuni Pozitie_strat[npart][n_max_vec_exch];

Moment M;

int nr_vecini[npart];
int nr_vecini_exch[npart];

double mhl;
double y_SW[nvar];
double y_in[nvar];
double y[nvar];
double amplitudini[npart];
///////////////////////////////////////////////////////
void init_particule(double th, double ph, double th_ea, double ph_ea) {
	ofstream fisier("E:\\Stoleriu\\C\\special\\3d\\res\\2021\\DSWM\\particule.dat");

	int j, counter;
	double dist = 2.0, dia;
	double V_inf;

	for (int i = 0; i < npart; i++) {
		P[i].index = i;
		P[i].R = 9.0 * R_mod;// pow(3.0*P[i].V / 4.0 / Pi, 1.0 / 3.0);
		P[i].V = 4.0 * Pi / 3.0 * P[i].R * P[i].R * P[i].R;// (1.0 + i * 0.0)*V_mod;// abs(V_mod*normal_rand(1.0, 0.05));
	}

	dia = 2.0 * P[0].R;
	V_inf = 4.0 * Pi / 3.0 * pow(dist * dia / 2.0, 3.0);

	for (int i = 0; i < npart; i++) {
		P[i].strat = 1;

		P[i].L_j = 5.0e-7 / V_inf;

		P[i].X = (i % 11) * dia * dist;// *uniform_rand(-1.001e-2, 1.001e-2);
		P[i].Y = (i / 11) * dia * dist;
		P[i].Z = 0.0 * dia * dist;// *uniform_rand(-1.001e-2, 1.001e-2);

		P[i].K1 = K1_mod;// *(1.5 + normal_rand(0.0001, 0.1));
		P[i].Ms = Ms_mod;// *(1.0 + normal_rand(0.0001, 0.1));

		P[i].theta_ea = gtr(th_ea);// +gtr(normal_rand(0.0001, 5.0));
		P[i].phi_ea = gtr(ph_ea);// +gtr(uniform_rand(0.0001, 359.0));

		P[i].T_relax = 1.0;// *(1.0 + normal_rand(0.0001, 0.1));

		P[i].theta_m = gtr(th);// +gtr(normal_rand(0.0, 5.0));
		P[i].phi_m = gtr(ph);// +gtr(uniform_rand(0.0, 359.0));
	}

	for (int i = 0; i < npart; i++) {
		P[i].K2 = 0.0;
		P[i].K3 = 0.0;

		P[i].m = P[i].Ms * P[i].V;
		P[i].Hk = 2.0 * P[i].K1 / miu0 / P[i].Ms;
		P[i].Hk2 = 4.0 * P[i].K2 / miu0 / P[i].Ms;
		P[i].T_Larmor_Hk = Pix2 / (gammap * abs(P[i].Hk)) * 1.0e+12;

		P[i].KV = P[i].K1 * P[i].V;

		P[i].Mx = P[i].Ms * sin(P[i].theta_m) * cos(P[i].phi_m);
		P[i].My = P[i].Ms * sin(P[i].theta_m) * sin(P[i].phi_m);
		P[i].Mz = P[i].Ms * cos(P[i].theta_m);

		P[i].eax = sin(P[i].theta_ea) * cos(P[i].phi_ea);
		P[i].eay = sin(P[i].theta_ea) * sin(P[i].phi_ea);
		P[i].eaz = cos(P[i].theta_ea);

		P[i].utheax = cos(P[i].theta_ea) * cos(P[i].phi_ea);
		P[i].utheay = cos(P[i].theta_ea) * sin(P[i].phi_ea);
		P[i].utheaz = -sin(P[i].theta_ea);

		P[i].upheax = -sin(P[i].phi_ea);
		P[i].upheay = cos(P[i].phi_ea);
		P[i].upheaz = 0.0;
	}

	Particula Data1, Data2;
	Interactiuni Pos_Coef;

	//initializare nr de vecini
	for (int i = 0; i < npart; i++) {
		nr_vecini[i] = nr_vecini_exch[i] = 0;
	}

	//alocarea vecinilor pentru interactiuni demagnetizante
	if (J_dip != 0.0) {
		for (int i = 0; i < npart; i++) {
			Data1 = P[i];
			for (int f = 0; f < npart; f++) {
				Data2 = P[f];
				if (i != f) {
					position_coeficients(Data1, Data2, &Pos_Coef);
					if (Pos_Coef.coef > P[i].L_j) {
						if (nr_vecini[i] + 1 > n_max_vec) {
							cout << "prea multi vecini de dipol" << endl;
							break;
						}
						else {
							nr_vecini[i]++;
							Pos_Coef.vecin = f;
							Pozitie[i][nr_vecini[i] - 1] = Pos_Coef;
						}
					}
				}
			}
		}
	}

	//alocarea vecinilor pentru exchange
	if (J_exch != 0.0) {
		for (int i = 0; i < npart; i++) {
			Data1 = P[i];
			for (int f = 0; f < npart; f++) {
				Data2 = P[f];
				if (i != f) {
					position_coeficients(Data1, Data2, &Pos_Coef);
					if (Pos_Coef.coef > prag_vecini_exch) {
						if (nr_vecini_exch[i] + 1 > n_max_vec_exch) {
							cout << "prea multi vecini de exchange" << endl;
							break;
						}
						else {
							nr_vecini_exch[i]++;
							Pos_Coef.vecin = f;
							Pozitie_strat[i][nr_vecini_exch[i] - 1] = Pos_Coef;
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < npart; i++)
		P_init[i] = P[i];

	fisier << "part" << '\t'
		<< "layer" << '\t'
		<< "X" << '\t'
		<< "Y" << '\t'
		<< "Z" << '\t'
		<< "K" << '\t'
		<< "th_EA" << '\t'
		<< "ph_EA" << '\t'
		<< "V" << '\t'
		<< "m" << '\t'
		<< "Ms" << '\t'
		<< "KV" << '\t'
		<< "vecini_dip" << '\t'
		<< "vecini_exch" << endl;

	for (int i = 0; i < npart; i++) {
		fisier << P[i].index << '\t'
			<< P[i].strat << '\t'
			<< P[i].X << '\t'
			<< P[i].Y << '\t'
			<< P[i].Z << '\t'
			<< P[i].K1 << '\t'
			<< rtg(P[i].theta_ea) << '\t'
			<< rtg(P[i].phi_ea) << '\t'
			<< P[i].V << '\t'
			<< P[i].m << '\t'
			<< P[i].Ms << '\t'
			<< P[i].V * P[i].K1 << '\t'
			<< nr_vecini[i] << '\t'
			<< nr_vecini_exch[i] << endl;
	}

	cout << "part" << " || "
		<< "layer" << " || "
		<< "X" << " || "
		<< "Y" << " || "
		<< "Z" << " || "
		<< "K" << " || "
		<< "th_EA" << " || "
		<< "ph_EA" << " || "
		<< "V" << " || "
		<< "m" << " || "
		<< "Ms" << " || "
		<< "KV" << " || "
		<< "vecini_dip" << " || "
		<< "vecini_exch" << endl;
	for (int i = 0; i < npart; i++) {
		cout << P[i].index << "  ||  "
			<< P[i].strat << "  ||  "
			<< P[i].X / dia << "  ||  "
			<< P[i].Y / dia << "  ||  "
			<< P[i].Z / dia << "  ||  "
			<< P[i].K1 / K1_mod << "  ||  "
			<< rtg(P[i].theta_ea) << "  ||  "
			<< rtg(P[i].phi_ea) << "  ||  "
			<< P[i].V << "  ||  "
			<< P[i].m << "  ||  "
			<< P[i].Ms / Ms_mod << "  ||  "
			<< P[i].KV << "  ||  "
			<< nr_vecini[i] << "  ||  "
			<< nr_vecini_exch[i] << "  ||  " << endl;
	}

	fisier.close();
	//cout << "////////S-AU GENERAT PARTICULELE/////" << endl;
}
///////////////////////////////////////////////////////////////////
void reset_particule() {
	for (int i = 0; i < npart; i++)
		P[i] = P_init[i];
}
/////////////////////////////////////////////////////////////
void Add_interactions(int part, Camp *H, double *y) {
	int kk;
	double mp_x, mp_y, mp_z;
	for (int j = 0; j < nr_vecini[part]; j++) {
		kk = 2 * Pozitie[part][j].vecin;
		mp_x = P[kk / 2].m * sin(y[kk + 0]) * cos(y[kk + 1]);
		mp_y = P[kk / 2].m * sin(y[kk + 0]) * sin(y[kk + 1]);
		mp_z = P[kk / 2].m * cos(y[kk + 0]);

		H->Hx += J_dip * Pozitie[part][j].coef * (Pozitie[part][j].xx * mp_x + Pozitie[part][j].xy * mp_y + Pozitie[part][j].xz * mp_z);
		H->Hy += J_dip * Pozitie[part][j].coef * (Pozitie[part][j].yx * mp_x + Pozitie[part][j].yy * mp_y + Pozitie[part][j].yz * mp_z);
		H->Hz += J_dip * Pozitie[part][j].coef * (Pozitie[part][j].zx * mp_x + Pozitie[part][j].zy * mp_y + Pozitie[part][j].zz * mp_z);
	}
	for (int j = 0; j < nr_vecini_exch[part]; j++) {
		kk = 2 * Pozitie_strat[part][j].vecin;
		mp_x = P[kk / 2].m * sin(y[kk + 0]) * cos(y[kk + 1]);
		mp_y = P[kk / 2].m * sin(y[kk + 0]) * sin(y[kk + 1]);
		mp_z = P[kk / 2].m * cos(y[kk + 0]);
		H->Hx += J_exch * mp_x;
		H->Hy += J_exch * mp_y;
		H->Hz += J_exch * mp_z;
	}
	H->Hamp = sqrt(H->Hx * H->Hx + H->Hy * H->Hy + H->Hz * H->Hz);
	H->theta_h = SafeAcos(H->Hz / H->Hamp);
	H->phi_h = atan2(H->Hy, H->Hx);
}
///////////////////////////////////////////////////////////////////
void build_magnetization(double *x, int lim1, int lim2, Moment *M) {
	M->Mx = 0.0;
	M->My = 0.0;
	M->Mz = 0.0;
	M->Mamp = 0.0;
	M->theta_M = 0.0;
	M->phi_M = 0.0;

	for (int i = lim1; i < lim2; i++) {
		P[i].theta_m = x[2 * i + 0];
		P[i].phi_m = x[2 * i + 1];

		P[i].Mx = P[i].Ms * sin(P[i].theta_m) * cos(P[i].phi_m);
		P[i].My = P[i].Ms * sin(P[i].theta_m) * sin(P[i].phi_m);
		P[i].Mz = P[i].Ms * cos(P[i].theta_m);

		M->Mx += P[i].Mx;
		M->My += P[i].My;
		M->Mz += P[i].Mz;

		M->Mamp += P[i].Ms;
	}

	M->theta_M = SafeAcos(M->Mz / M->Mamp);
	M->phi_M = atan2(M->My / M->Mamp, M->Mx / M->Mamp);
	M->Mx = sin(M->theta_M) * cos(M->phi_M);
	M->My = sin(M->theta_M) * sin(M->phi_M);
	M->Mz = cos(M->theta_M);
}
///////////////////////////////////////////////////////////////////
double SW_fcn(double x, void *params) {
	struct SW_function_params *p = (struct SW_function_params *)params;

	double Hk_fcn = p->Hk_fcn;
	double H_fcn = p->H_fcn;
	double theta0_fcn = p->theta0_fcn;

	return Hk_fcn * cos(x) * sin(x) + H_fcn * sin(x - theta0_fcn);
}
///////////////////////////////////////////////////////////////////
void SW(double th_h, double ph_h, double a, double *sol) {
	double st, sp, ct, cp;
	double mx, my, mz;
	double hx, hy, hz;
	double loco_angle_H, loco_angle_M, loco_old_angle_M, altr_angle_M;
	double its_a_sin, its_a_cos, aux1, aux2, g, Hc;
	double h_w, hx_w, hy_w, d;
	double temp_Hea, temp_Hha, temp_angle_H, temp_ea, temp_H;

	double Sol_bune[2], sol_bune[2], x_lo, x_hi, r;

	double *sol_j = new double[nvar];
	for (int part = 0; part < npart; part++) {
		sol_j[2 * part + 0] = sol[2 * part + 0];
		sol_j[2 * part + 1] = sol[2 * part + 1];
	}

	for (int part = 0; part < npart; part++) {
		camp(th_h, ph_h, a, &H);
		Add_interactions(part, &H, sol_j);
		amplitudini[part] = H.Hamp;

		st = sin(sol[2 * part + 0]);
		sp = sin(sol[2 * part + 1]);
		ct = cos(sol[2 * part + 0]);
		cp = cos(sol[2 * part + 1]);

		mx = st * cp;
		my = st * sp;
		mz = ct;

		hx = sin(H.theta_h) * cos(H.phi_h);
		hy = sin(H.theta_h) * sin(H.phi_h);
		hz = cos(H.theta_h);

		loco_angle_H = SafeAcos(P[part].eax * hx + P[part].eay * hy + P[part].eaz * hz);
		loco_old_angle_M = SafeAcos(P[part].eax * mx + P[part].eay * my + P[part].eaz * mz);

		its_a_sin = sin(loco_angle_H);
		its_a_cos = cos(loco_angle_H);
		aux1 = pow(its_a_sin * its_a_sin, (1.0 / 3.0));
		aux2 = pow(its_a_cos * its_a_cos, (1.0 / 3.0));
		g = pow(aux1 + aux2, -1.5);
		Hc = P[part].Hk * g;

		h_w = H.Hamp / abs(P[part].Hk);
		hy_w = h_w * its_a_cos;
		hx_w = h_w * its_a_sin;
		d = 1 - h_w * h_w;

		doublec e(0.0, 0.0), fp(0.0, 0.0), fm(0.0, 0.0), mp(0.0, 0.0), mm(0.0, 0.0), mpx(0.0, 0.0), mmx(0.0, 0.0), t1(0.0, 0.0), t2(0.0, 0.0);

		e = d * cos(acos(doublec(54.0 * hx_w * hx_w * hy_w * hy_w / d / d / d - 1.0)) / 3.0);
		fp = sqrt(9.0 * hy_w * hy_w + 6.0 * d + 6.0 * e);
		fm = -fp;
		mp = (fp + sqrt(2.0 * fp * fp - 18.0 * e + 54.0 * hy_w * (1.0 + hx_w * hx_w) / fp)) / 6.0 - hy_w / 2.0;
		mm = (fm - sqrt(2.0 * fm * fm - 18.0 * e + 54.0 * hy_w * (1.0 + hx_w * hx_w) / fm)) / 6.0 - hy_w / 2.0;
		mpx = sqrt(1.0 - mp * mp);
		mmx = sqrt(1.0 - mm * mm);
		t1 = atan(mpx / mp);
		t2 = atan(mmx / mm);
		sol_bune[0] = real(t1);
		sol_bune[1] = Pi + real(t2);

		if (H.Hamp == 0.0) {
			loco_angle_M = (loco_old_angle_M < Pis2) ? 0.0 : Pi;
		}
		if (H.Hamp >= Hc) {
			loco_angle_M = (loco_angle_H < Pis2) ? sol_bune[0] : sol_bune[1];
		}
		if (H.Hamp < Hc) {
			loco_angle_M = (fabs(loco_old_angle_M - sol_bune[0]) > fabs(loco_old_angle_M - sol_bune[1])) ? sol_bune[1] : sol_bune[0];
			altr_angle_M = (fabs(loco_old_angle_M - sol_bune[0]) > fabs(loco_old_angle_M - sol_bune[1])) ? sol_bune[0] : sol_bune[1];
		}

		temp_ea = (loco_angle_M <= Pi) ? sin(loco_angle_H - loco_angle_M) / sin(loco_angle_H) : sin(loco_angle_H - (Pix2 - loco_angle_M)) / sin(loco_angle_H);
		temp_H = (loco_angle_M <= Pi) ? sin(loco_angle_M) / sin(loco_angle_H) : sin(Pix2 - loco_angle_M) / sin(loco_angle_H);

		mx = (temp_ea * P[part].eax + temp_H * hx);
		my = (temp_ea * P[part].eay + temp_H * hy);
		mz = (temp_ea * P[part].eaz + temp_H * hz);

		sol[2 * part + 0] = SafeAcos(mz);
		sol[2 * part + 1] = atan2(my, mx);
	}
	delete[] sol_j;
	/*
		gsl_function F;
	int n_sols, status, iter, max_iter = 203;

	for (int part = 0; part < npart; part++) {
		iter = 0;
		double *sols, *limits;
		const gsl_root_fsolver_type *T;
		gsl_root_fsolver *s;
		struct SW_function_params params = { P[part].Hk, H_SW[part].Hamp, loco_angle_H[part] };
		F.function = &SW_fcn;
		F.params = &params;
		T = gsl_root_fsolver_bisection;

		if (H_SW[part].Hamp == 0.0) {
			sols = (double *)calloc(1, sizeof(double));
			limits = (double *)calloc(1, sizeof(double));
			loco_angle_M[part] = (loco_old_angle_M[part] < Pis2) ? 0.0 : Pi;
		}
		else {
			temp_Hea = H_SW[part].Hamp * cos(loco_angle_H[part]);
			temp_Hha = H_SW[part].Hamp * sin(loco_angle_H[part]);
			if ((temp_Hea >= 0.0) && (temp_Hha >= 0.0))
				temp_angle_H = Pi + atan(-pow(tan(loco_angle_H[part])*tan(loco_angle_H[part]), (1.0 / 6.0)));	//C1
			if ((temp_Hea >= 0.0) && (temp_Hha < 0.0))
				temp_angle_H = Pi + atan(pow(tan(loco_angle_H[part])*tan(loco_angle_H[part]), (1.0 / 6.0)));	//C4
			if ((temp_Hea < 0.0) && (temp_Hha < 0.0))
				temp_angle_H = Pix2 + atan(-pow(tan(loco_angle_H[part])*tan(loco_angle_H[part]), (1.0 / 6.0)));	//C3
			if ((temp_Hea < 0.0) && (temp_Hha >= 0.0))
				temp_angle_H = atan(pow(tan(loco_angle_H[part])*tan(loco_angle_H[part]), (1.0 / 6.0)));	//C2

			n_sols = (H_SW[part].Hamp >= Hc[part]) ? 2 : 4;
			sols = (double *)calloc(n_sols, sizeof(double));
			limits = (double *)calloc(2 * n_sols, sizeof(double));

			if (H_SW[part].Hamp >= Hc[part]) {
				if ((temp_Hea >= 0.0) && (temp_Hha >= 0.0)) {
					limits[0] = 0.0; limits[1] = Pis2; limits[2] = Pi; limits[3] = 3.0*Pis2;	//C1
				}
				if ((temp_Hea >= 0.0) && (temp_Hha < 0.0)) {
					limits[0] = Pis2; limits[1] = Pi; limits[2] = 3.0*Pis2; limits[3] = Pix2;	//C4
				}
				if ((temp_Hea < 0.0) && (temp_Hha < 0.0)) {
					limits[0] = 0.0; limits[1] = Pis2; limits[2] = Pi; limits[3] = 3.0*Pis2;	//C3
				}
				if ((temp_Hea < 0.0) && (temp_Hha >= 0.0)) {
					limits[0] = Pis2; limits[1] = Pi; limits[2] = 3.0*Pis2; limits[3] = Pix2;	//C2
				}
			}
			else {
				if ((temp_Hea >= 0.0) && (temp_Hha >= 0.0)) {
					limits[0] = 0.0; limits[1] = Pis2; limits[2] = Pi; limits[3] = 3.0*Pis2; limits[4] = Pis2; limits[5] = temp_angle_H; limits[6] = temp_angle_H; limits[7] = Pi;
				}
				if ((temp_Hea >= 0.0) && (temp_Hha < 0.0)) {
					limits[0] = Pis2; limits[1] = Pi; limits[2] = 3.0*Pis2; limits[3] = Pix2; limits[4] = Pi; limits[5] = temp_angle_H; limits[6] = temp_angle_H; limits[7] = 3.0*Pis2;
				}
				if ((temp_Hea < 0.0) && (temp_Hha < 0.0)) {
					limits[0] = 0.0; limits[1] = Pis2; limits[2] = Pi; limits[3] = 3.0*Pis2; limits[4] = 3.0*Pis2; limits[5] = temp_angle_H; limits[6] = temp_angle_H; limits[7] = Pix2;
				}
				if ((temp_Hea < 0.0) && (temp_Hha >= 0.0)) {
					limits[0] = Pis2; limits[1] = Pi; limits[2] = 3.0*Pis2; limits[3] = Pix2; limits[4] = 0.0; limits[5] = temp_angle_H; limits[6] = temp_angle_H; limits[7] = Pis2;
				}
			}

			for (int scontor = 0; scontor < n_sols; scontor++) {
				s = gsl_root_fsolver_alloc(T);
				x_lo = limits[2 * scontor + 0];
				x_hi = limits[2 * scontor + 1];
				gsl_root_fsolver_set(s, &F, x_lo, x_hi);
				do {
					iter++;
					status = gsl_root_fsolver_iterate(s);
					r = gsl_root_fsolver_root(s);
					x_lo = gsl_root_fsolver_x_lower(s);
					x_hi = gsl_root_fsolver_x_upper(s);
					status = gsl_root_test_interval(x_lo, x_hi, 0, 1.0e-7);
				} while (status == GSL_CONTINUE && iter < max_iter);
				sols[scontor] = r;
				gsl_root_fsolver_free(s);
			}

			if (H_SW[part].Hamp >= Hc[part]) {
				for (int scontor = 0; scontor < n_sols; scontor++)
					if ((P[part].Hk * cos(2.0*sols[scontor]) + H_SW[part].Hamp * cos(sols[scontor] - loco_angle_H[part])) > 0)
						loco_angle_M[part] = sols[scontor];
			}
			else {
				int contor_sol_bune = 0;
				for (int scontor = 0; scontor < n_sols; scontor++)
					if ((P[part].Hk * cos(2.0*sols[scontor]) + H_SW[part].Hamp * cos(sols[scontor] - loco_angle_H[part])) > 0) {
						Sol_bune[contor_sol_bune] = sols[scontor];
						contor_sol_bune++;
					}
				loco_angle_M[part] = (fabs(loco_old_angle_M[part] - Sol_bune[0]) > fabs(loco_old_angle_M[part] - Sol_bune[1])) ? Sol_bune[1] : Sol_bune[0];
			}
		}
	}*/
}
///////////////////////////////////////////////////////////////////

void dynamic_SW(double t, double *sol_old, double *sol_target, double *sol_new) 
{
	double raport;
	double mx0, mx1, my0, my1, mz0, mz1, mxf, myf, mzf;
	double mx0r, my0r, mz0r, mx1r, my1r, mz1r;
	double t_max, loco_angle_with_ea_old, loco_angle_with_ea_target, angle_target, angle_new;
	double th_m0r, th_m1r, ph_m0r, ph_m1r;
	double s0, c0, s1, c1;
	double loco_angle_old_uthea, loco_angle_target_uthea;
	double torqx, torqy, torqz;

	for (int i = 0; i < npart; i++) {
		mx0 = sin(sol_old[2 * i + 0]) * cos(sol_old[2 * i + 1]);
		my0 = sin(sol_old[2 * i + 0]) * sin(sol_old[2 * i + 1]);
		mz0 = cos(sol_old[2 * i + 0]);

		mx1 = sin(sol_target[2 * i + 0]) * cos(sol_target[2 * i + 1]);
		my1 = sin(sol_target[2 * i + 0]) * sin(sol_target[2 * i + 1]);
		mz1 = cos(sol_target[2 * i + 0]);

		loco_angle_with_ea_old = SafeAcos(P[i].eax * mx0 + P[i].eay * my0 + P[i].eaz * mz0);
		angle_target = SafeAcos(mx0 * mx1 + my0 * my1 + mz0 * mz1);
		loco_angle_with_ea_target = SafeAcos(P[i].eax * mx1 + P[i].eay * my1 + P[i].eaz * mz1);
		loco_angle_old_uthea = SafeAcos(mx0 * P[i].utheax + my0 * P[i].utheay + mz0 * P[i].utheaz);
		loco_angle_target_uthea = SafeAcos(mx1 * P[i].utheax + my1 * P[i].utheay + mz1 * P[i].utheaz);

		if ((loco_angle_target_uthea<Pis2 && loco_angle_old_uthea>Pis2) || (loco_angle_target_uthea > Pis2 && loco_angle_old_uthea < Pis2))
			loco_angle_with_ea_old = Pix2 - loco_angle_with_ea_old;

		// se roteste pozitia veche a lui H astfel incat acesta sa fie dupa Oz si se afla noul M
		double R1[3][3], R2[3][3];
		////////////////////////////////////////////////////////////////////////// 
		////////////////////////////////////////////////////////////////////////// CE-I ASTA???

		s0 = sin(sol_target[2 * i + 0]);
		s1 = sin(sol_target[2 * i + 1]);
		c0 = cos(sol_target[2 * i + 0]);
		c1 = cos(sol_target[2 * i + 1]);

		R1[0][0] = c0 * c1;  R1[0][1] = c0 * s1;  R1[0][2] = -s0;
		R1[1][0] = -s1;       R1[1][1] = c1;       R1[1][2] = 0.0;
		R1[2][0] = s0 * c1;  R1[2][1] = s0 * s1;  R1[2][2] = c0;

		R2[0][0] = R1[0][0]; R2[0][1] = R1[1][0]; R2[0][2] = R1[2][0];
		R2[1][0] = R1[0][1]; R2[1][1] = R1[1][1]; R2[1][2] = R1[2][1];
		R2[2][0] = R1[0][2]; R2[2][1] = R1[1][2]; R2[2][2] = R1[2][2];

		mx0r = R1[0][0] * mx0 + R1[0][1] * my0 + R1[0][2] * mz0;
		my0r = R1[1][0] * mx0 + R1[1][1] * my0 + R1[1][2] * mz0;
		mz0r = R1[2][0] * mx0 + R1[2][1] * my0 + R1[2][2] * mz0;

		th_m0r = SafeAcos(mz0r);
		ph_m0r = atan2(my0r, mx0r);

		// se modifica th_m0r si ph_m0r conform regulilor
		raport = amplitudini[i] / P[i].Hk;
		raport = (raport < 2.0) ? 2.0 : raport;
		t_max = P[i].T_Larmor_Hk * timp_max(raport, loco_angle_with_ea_old, loco_angle_with_ea_target) / P[i].T_relax;

		if (angle_target < Pis2) {
			angle_new = 4.0 * angle_target * (1.0 + fabs(cos(angle_target)));
		}
		else {
			angle_new = -4.0 * (angle_target - Pi) * (1.0 + fabs(cos(angle_target)));
		}

		th_m1r = (t > t_max) ? 0.0 : th_m0r - angle_new * t / t_max;
		ph_m1r = ph_m0r - Pix2 * t * P[i].T_relax * uniform_rand(1.0, 1.13) / P[i].T_Larmor_Hk;

		// se roteste totul inapoi rezultand noua solutie
		mx1r = sin(th_m1r) * cos(ph_m1r);
		my1r = sin(th_m1r) * sin(ph_m1r);
		mz1r = cos(th_m1r);

		mxf = R2[0][0] * mx1r + R2[0][1] * my1r + R2[0][2] * mz1r;
		myf = R2[1][0] * mx1r + R2[1][1] * my1r + R2[1][2] * mz1r;
		mzf = R2[2][0] * mx1r + R2[2][1] * my1r + R2[2][2] * mz1r;

		torqx = myf * mz1 - mzf * my1;
		torqy = mzf * mx1 - mxf * mz1;
		torqz = mxf * my1 - myf * mx1;

		torq_mod = sqrt(torqx * torqx + torqy * torqy + torqz * torqz);

		sol_new[2 * i + 0] = SafeAcos(mzf);
		sol_new[2 * i + 1] = atan2(myf, mxf);
	}
}
///////////////////////////////////////////////////////////////////

void mhl_sw(double th_h, double ph_h, double t_stop, double T, double pas_SW) 
{
	cout << "|||| SE CALCULEAZA MHL CU SW |||" << endl;
	ofstream fisier("E:\\Stoleriu\\C\\special\\3d\\res\\2021\\DSWM\\mhl_sw.dat");
	double h = 0, t = 0;

	for (int i = 0; i < (int)(sir_camp.size()); i++) {
		h = sir_camp[i];
		t = 0.0;
		torq_mod = 1.0;

		while (cond_sw) 
		{
			for (int part = 0; part < npart; part++) 
			{
				y_in[2 * part + 0] = P[part].theta_m;
				y_in[2 * part + 1] = P[part].phi_m;
				y[2 * part + 0] = y_in[2 * part + 0];
				y[2 * part + 1] = y_in[2 * part + 1];
			}

			SW(th_h, ph_h, h, y);
			
			for (int part = 0; part < npart; part++) 
			{
				y_SW[2 * part + 0] = y[2 * part + 0];
				y_SW[2 * part + 1] = y[2 * part + 1];
				y[2 * part + 0] = y_in[2 * part + 0];
				y[2 * part + 1] = y_in[2 * part + 1];
			}

			dynamic_SW(pas_SW, y_in, y_SW, y);
			
			for (int part = 0; part < npart; part++) 
			{
				P[part].theta_m = y[2 * part + 0];
				P[part].phi_m = y[2 * part + 1];
			}
			t += pas_SW;
		}

		camp(th_h, ph_h, h, &H);
		build_magnetization(y, 0, npart, &M);
		mhl = (M.Mx * H.Hx + M.My * H.Hy + M.Mz * H.Hz) / h;
		cout << i << '\t' << h / 1 << '\t' << mhl << endl;
		fisier << h / 1 << '\t' << mhl << endl;
	}
	fisier.close();
}
///////////////////////////////////////////////////
int main()
{
	srand((unsigned)time(0));

	double mu_H = 10.001;//Frecventa campului extern in MHz
	double mu_punct = mu_H * 405;//Frecventa de inregistrare a punctelor de pe MHL in MHz
	double T_H = 1.0 / mu_H * 1.0e+6;//ps
	double t_stop = 1.0 / mu_punct * 1.0e+6;//ps

	cout << "t_stop = " << t_stop << " picosecunde" << endl;
	cout << "T_H = " << T_H << " picosecunde" << endl;

	init_particule(0.0015, 1.5e-2, 0.12, 0.12);
	camp_sir(0, 600000, t_stop, T_H, &sir_camp);
	mhl_sw(0.011, 0.012, t_stop, T_H, t_stop / 250.0);

	return 0;
}