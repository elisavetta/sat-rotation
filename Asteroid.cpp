#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <array>
#include "Integrator.h"

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using namespace std;
typedef std::array<double, 5> OutArray;

/*
const std::size_t Dimension = 8; // solving 8 equations 
using Array = std::array<double, Dimension>;
typedef std::array<double, 8> OutArray;
*/

const double PI = 3.14159265358979323846;
const double twoPI = 2 * PI;
const double deg2rad = PI / 180;
const double years2sec = 365.25 * 24 * 3600;

int sign(double a)
{
	return ((a > 0) - (a < 0));
}

Vector3d rotateByQ(const Vector4d& q, const Vector3d& r)
{
	Vector3d x;
	x(0) = (q(0) * q(0) + q(1) * q(1) - q(2) * q(2) - q(3) * q(3)) * r(0) + 2 * (q(1) * q(2) + q(0) * q(3)) * r(1) + 2 * (q(1) * q(3) - q(0) * q(2)) * r(2);
	x(1) = 2 * (q(1) * q(2) - q(0) * q(3)) * r(0) + (q(0) * q(0) + q(2) * q(2) - q(1) * q(1) - q(3) * q(3)) * r(1) + 2 * (q(2) * q(3) + q(0) * q(1)) * r(2);
	x(2) = 2 * (q(1) * q(3) + q(0) * q(2)) * r(0) + 2 * (q(2) * q(3) - q(0) * q(1)) * r(1) + (q(0) * q(0) + q(3) * q(3) - q(1) * q(1) - q(2) * q(2)) * r(2);
	return x;
}

Vector4d quatproduct(const Vector4d& q, const Vector4d& m) {
	Vector4d x; // function is used for finding initial quaternion
	x(0) = q(0) * m(0) - q(1) * m(1) - q(2) * m(2) - q(3) * m(3);
	x(1) = q(0) * m(1) + q(1) * m(0) + q(2) * m(3) - q(3) * m(2);
	x(2) = q(0) * m(2) + q(2) * m(0) + q(3) * m(1) - q(1) * m(3);
	x(3) = q(0) * m(3) + q(3) * m(0) + q(1) * m(2) - q(2) * m(1);
	return x;
}

struct orbparam
{
	double mu, a_major, eccentricity, w0, w02;
};

struct satparam
{
	satparam()
	{
		itm_B(0, 0) = 300; itm_B(1, 1) = 200; itm_B(2, 2) = 100;  // inertia tensor
		itm_B(0, 1) = itm_B(0, 2) = itm_B(1, 0) = itm_B(1, 2) = itm_B(2, 0) = itm_B(2, 1) = 0;
	}
	Eigen::Matrix3d itm_B;  //inertia tensor matrix
};

void rhside(const satparam& sat, const orbparam& orb, const double t, const Array& x, Array* dx)
{
	Vector4d q = { x[0], x[1], x[2], x[3] };
	Vector3d omega = { x[4], x[5], x[6] };
	double nu = x[7]; //true anomaly

	Vector3d omega_dot;

	//recalculate OY and OZ unit vectors of the orbital reference frame wrt PG frame
	Vector3d eRorb_B, eRorb_PG;
	eRorb_PG = { cos(x[7]), sin(x[7]), 0};
	eRorb_B = rotateByQ(q, eRorb_PG);

	double ecc_param3 = (1 - orb.eccentricity * orb.eccentricity) * (1 - orb.eccentricity * orb.eccentricity) * (1 - orb.eccentricity * orb.eccentricity);
	double eccos = (1 + orb.eccentricity * cos(x[7]));

	//gravity gradient torque
	Vector3d Tgg;
	Tgg = 3 * (orb.w02) * eccos * eccos * eccos * eRorb_B.cross(sat.itm_B * eRorb_B)/ ecc_param3;

	//gyroscpoic term
	Vector3d wJw;
	wJw = omega.cross(sat.itm_B * omega);
	omega_dot = sat.itm_B.inverse() * (-wJw + Tgg);

	//true anomaly term
	double dnu = orb.w0 / (eccos *sqrt(ecc_param3)); 

	(*dx)[0] = 0.5 * (-q[1] * omega[0] - q[2] * omega[1] - q[3] * omega[2]);
	(*dx)[1] = 0.5 * (q[0] * omega[0] + q[2] * omega[2] - q[3] * omega[1]);
	(*dx)[2] = 0.5 * (q[0] * omega[1] - q[1] * omega[2] + q[3] * omega[0]);
	(*dx)[3] = 0.5 * (q[0] * omega[2] + q[1] * omega[1] - q[2] * omega[0]);

	(*dx)[4] = omega_dot[0];
	(*dx)[5] = omega_dot[1];
	(*dx)[6] = omega_dot[2];
	(*dx)[7] = dnu;

}

double YacobiIntegral(const orbparam& orb, const satparam& sat, double t, const Array& x)
{
	Vector4d q = { x[0], x[1], x[2], x[3] };
	Vector3d omega = { x[4], x[5], x[6] };
	double nu = x[7];

	//get OX and OZ unit vectors of the orbital reference frame wrt PG frame
	Vector3d eRorb_B, eRorb_PG; //  OX unit vectors of the orbital reference frame wrt PG frame
	eRorb_PG = { cos(x[7]), sin(x[7]), 0 };
	eRorb_B = rotateByQ(q, eRorb_PG);

	Vector3d eZ_PG, eZ_B;
	eZ_PG = { 0, 0, 1 };
	eZ_B = rotateByQ(q, eZ_PG);

	Vector3d wt = omega - orb.w0 * eZ_B; 
	
	double Tr, Vg, T0;

	Tr = 0.5 * wt.transpose() * (sat.itm_B * wt);
	Vg = 1.5 * orb.w02 * eRorb_B.transpose() * (sat.itm_B * eRorb_B);
	T0 = 0.5 * orb.w02 * eZ_B.transpose() * (sat.itm_B * eRorb_B);

	return Tr + Vg - T0;
}

void output(const satparam& sat, const orbparam& orb, const double t, const Array& x, OutArray* out)
{
	Vector4d q = { x[0], x[1], x[2], x[3] };
	Vector3d omega = { x[4], x[5], x[6] };

	Vector3d K_B, eK_B;
	K_B = sat.itm_B * omega;
	eK_B = K_B / sqrt(K_B(0) * K_B(0) + K_B(1) * K_B(1) + K_B(2) * K_B(2));

	Vector3d eXpg_B, eYpg_B, eZpg_B;
	eXpg_B = rotateByQ(q, { 1, 0, 0 });
	eYpg_B = rotateByQ(q, { 0, 1, 0 });
	eZpg_B = rotateByQ(q, { 0, 0, 1 });

	double rho = acos(eK_B.dot(eZpg_B));
	double sin_rho = sin(rho);
	double tetta = acos(eK_B(2));
	double sigma_s, cos_sigma_s, sin_sigma_s;
	
	if (abs(sin_rho) <= 1e-5)
	{
		sigma_s = 0;
	}
	else
	{
		cos_sigma_s = (eK_B.dot(eXpg_B)) / sin_rho;
		sin_sigma_s = (eK_B.dot(eYpg_B)) / sin_rho;

		if (abs(cos_sigma_s) > 1)
			sigma_s = 0.5 * PI * (1 - sign(cos_sigma_s));
		else
			if (sin_sigma_s < 0)
				sigma_s = PI - acos(cos_sigma_s);
			else
				sigma_s = acos(cos_sigma_s);
	}

	Vector3d el1_pg, el1_B, el2_pg, el2_B;
	//el1_SO = { cos(rho) * cos(sigma_s), cos(rho) * sin(sigma_s), -sin(rho) }; 
	//el1_B = rotateByQ(q, el1_pg); //wrong value
	el2_pg = { -sin(sigma_s), cos(sigma_s), 0 };
	el2_B = rotateByQ(q, el2_pg);
	el1_B = el2_B.cross(eK_B);

	double psi;
	if (abs(sin(tetta)) <= 1e-5)
		psi = acos(el1_B(0));
	else
		psi = atan(-el1_B(2) / eK_B(2));

	(*out)[0] = t / 86400;
	(*out)[1] = sigma_s / deg2rad;
	(*out)[2] = rho / deg2rad;
	(*out)[3] = psi / deg2rad;
	(*out)[4] = tetta / deg2rad;
	/*
	(*out)[2] = el1_B(0);
	(*out)[3] = el1_B(1);
	(*out)[4] = el1_B(2);
	*/
}

int main()
{
	satparam sat;
	orbparam orb;
	Integrator integrator;
	Array x;
	OutArray x_out;

	orb.mu = 3.986e14;     //the Earth gravity constant
	orb.a_major = 1340000; //semi-major axis
	orb.eccentricity = 0;
	orb.w0 = sqrt(orb.mu)/sqrt(orb.a_major* orb.a_major* orb.a_major);
	orb.w02 = orb.w0 * orb.w0;

	//Initial conditions
	double sigma_0, ro_0, psi_0, fi_0, tetta_0;
	sigma_0 = 30*deg2rad;
	ro_0 = 120 * deg2rad;
	psi_0 = 0;
	tetta_0 = 0;
	fi_0 = 0;
	double mod_omega_0 = 1; //module of omega_0
	double nu_0 = 0; //initial value of true anomaly

	//initial value of omega
	Vector3d omega_0 = { mod_omega_0*sin(tetta_0)*sin(fi_0), mod_omega_0*sin(tetta_0)*cos(fi_0), mod_omega_0*cos(tetta_0) };

	//calculation of the initial quaternion value
	Vector4d q_1, q_2, q_3, q_4, q_5, q_6, q_7, q_8, q_0;
	q_1 = { cos(sigma_0 / 2), 0, 0, sin(sigma_0 / 2) };
	q_2 = { cos(ro_0 / 2), 0, -sin(ro_0 / 2), 0 };
	q_3 = quatproduct(q_1, q_2);
	q_4 = { cos(psi_0 / 2), 0, 0, sin(psi_0 / 2) };
	q_5 = quatproduct(q_3, q_4);
	q_6 = { cos(tetta_0 / 2), sin(tetta_0 / 2), 0, 0};
	q_7 = quatproduct(q_5, q_6);
	q_8 = { cos(fi_0 / 2), 0, 0, sin(fi_0 / 2) };
	q_0 = quatproduct(q_7, q_8); // initial quaternion value

	//initial condition vector
	x[0] = q_0[0];
	x[1] = q_0[1];
	x[2] = q_0[2]; 
	x[3] = q_0[3];
	x[4] = omega_0[0];
	x[5] = omega_0[1];
	x[6] = omega_0[2];
	x[7] = nu_0;

	//integrator parameters
	double t = 0;
	double h = 1;
	double t_end = 0.4 * years2sec;
	double t_save = 100;
	int save_counter = 0;

	std::cout.precision(6);
	std::cout.setf(std::ios::fixed);
	std::ofstream  resultFile("object.csv");
	
	//cout << t << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << endl;
	//cout << t << "," << YacobiIntegral(orb, sat, t, x) << endl;

	while (t < t_end)
	{
		double q_norm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
		for (int i = 0; i < 4; i++) x[i] = x[i] / q_norm;
		if (t >= save_counter * t_save)
		{
			output(sat, orb, t, x, &x_out);
			resultFile << setprecision(8) << x_out[0] << "," << x_out[1] << "," << x_out[2] << "," << x_out[3] << "," << x_out[4] << "," << YacobiIntegral(orb, sat, t, x) << endl;
			//cout << t << "," << YacobiIntegral(orb, sat, t, x) << endl;
			cout << x_out[0] << " " << x_out[1] << " " << x_out[2] << " " << x_out[3] << " " << x_out[4] << " " << YacobiIntegral(orb, sat, t, x) << endl;
			save_counter += 1;
		}
		integrator.rk78(&t, &x, &h, 1e-9, 1e-3, 1e+3, [&sat, &orb](const double t, const Array& x, Array* dx) { rhside(sat, orb, t, x, dx); });
		//cout << t << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << endl;

	}

	double q_norm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
	for (int i = 0; i < 4; i++) x[i] = x[i] / q_norm;

	output(sat, orb, t, x, &x_out);
	resultFile << setprecision(8) << x_out[0] << "," << x_out[1] << "," << x_out[2] << "," << x_out[3] << "," << x_out[4] << "," << YacobiIntegral(orb, sat, t, x) << endl;
	resultFile.close();

	return 0;
}