#pragma once
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define SGN(a) (((a)<0) ? -1 : 1)

const std::size_t Dimension = 8; // solving 8 equations 
using Array = std::array<double, Dimension>;

using Eigen::Matrix3d;
using Eigen::Vector3d;
using namespace std;


//integrator parameters
struct integParams {
	int n = 8;         // number of equations
	double h = 1;               // initial step, in seconds
	double hMin = 1e-3;         // min step, in seconds
	double hMax = 1e+3;         // max step, in seconds
	double tolerance = 1e-9;    // integrator tolerance
	double t0 = 0;              // start time 
	double t_end = 10;           // end time, in years
	double t_save = 100;           // save output every t_save, in seconds
};

static const std::array<double, 13> alfa = {
    0.e0,      2.e0 / 27.e0,       1.e0 / 9.e0,        1.e0 / 6.e0,
    5.e0 / 12.e0,           0.5e0,       5.e0 / 6.e0,        1.e0 / 6.e0,
    2.e0 / 3.e0,       1.e0 / 3.e0,            1.e0,             0.e0,
    1.e0 };

static const std::array<double, 79> beta = {
    0.e0,      2.e0 / 27.e0,      1.e0 / 36.e0,       1.e0 / 12.e0,
    1.e0 / 24.e0,            0.e0,       1.e0 / 8.e0,       5.e0 / 12.e0,
    0.e0,    -25.e0 / 16.e0,     25.e0 / 16.e0,            .5e-1,
    0.e0,            0.e0,           .25e0,             .2e0,
    -25.e0 / 108.e0,            0.e0,            0.e0,    125.e0 / 108.e0,
    -65.e0 / 27.e0,    125.e0 / 54.e0,    31.e0 / 300.e0,             0.e0,
    0.e0,            0.e0,    61.e0 / 225.e0,       -2.e0 / 9.e0,
    13.e0 / 900.e0,            2.e0,            0.e0,             0.e0,
    -53.e0 / 6.e0,    704.e0 / 45.e0,    -107.e0 / 9.e0,      67.e0 / 90.e0,
    3.e0,   -91.e0 / 108.e0,            0.e0,             0.e0,
    23.e0 / 108.e0,  -976.e0 / 135.e0,    311.e0 / 54.e0,     -19.e0 / 60.e0,
    17.e0 / 6.e0,     -1.e0 / 12.e0, 2383.e0 / 4100.e0,             0.e0,
    0.e0,  -341.e0 / 164.e0, 4496.e0 / 1025.e0,    -301.e0 / 82.e0,
    2133.e0 / 4100.e0,     45.e0 / 82.e0,    45.e0 / 164.e0,      18.e0 / 41.e0,
    3.e0 / 205.e0,            0.e0,            0.e0,             0.e0,
    0.e0,     -6.e0 / 41.e0,    -3.e0 / 205.e0,      -3.e0 / 41.e0,
    3.e0 / 41.e0,      6.e0 / 41.e0,            0.e0, -1777.e0 / 4100.e0,
    0.e0,            0.e0,  -341.e0 / 164.e0,  4496.e0 / 1025.e0,
    -289.e0 / 82.e0, 2193.e0 / 4100.e0,     51.e0 / 82.e0,     33.e0 / 164.e0,
    12.e0 / 41.e0,            0.e0,            1.e0 };

static const std::array<double, 11> c7 = {
    41.e0 / 840.e0,            0.e0,            0.e0,             0.e0,
    0.e0,    34.e0 / 105.e0,      9.e0 / 35.e0,       9.e0 / 35.e0,
    9.e0 / 280.e0,     9.e0 / 280.e0,    41.e0 / 840.e0 };

static const std::array<double, 13> c8 = {
    0.e0,            0.e0,            0.e0,             0.e0,
    0.e0,    34.e0 / 105.e0,      9.e0 / 35.e0,       9.e0 / 35.e0,
    9.e0 / 280.e0,     9.e0 / 280.e0,            0.e0,     41.e0 / 840.e0,
    41.e0 / 840.e0 };


class Integrator
{
private:
    std::array<Array, 13> k;
    Array x7, x8, xpon, dx;

public:
    Integrator()
    {
        Array x;
        static_assert(x.size() >= 1, "");
    }

double rk78(double *at, Array* x, double *ah, double tol, double hmin, double hmax, const std::function<void(const double t, const Array& x, Array* dx)>& field)

	/*
	this routine performs one step of the integration procedure.
	the initial condition (at,x) is changed by a new one corresponding
	to the same orbit. the error is controlled by the threshold tol,
	and an estimate of the error produced in the actual step is returned
	as the value of the function.

	parameters:
	at:   time. input: time corresponding to the actual initial condition.
	output: new value corresponding to the new initial condition.
	x:    position. same remarks as at.
	ah:   time step (it can be modified by the routine according to the
	given threshold).
	tol:  threshold to control the integration error.
	hmin: minimun time step allowed.
	hmax: maximum time step allowed.
	n:    dimension of the system of odes.
	field: function that returns the value of the vector field.
	paramSet: set of parameters

	returned value: an estimate of the error produced in the actual step of
	integration.
	*/
{
	double tpon, tol1, err, nor, h1;
	
	do {
		
		//this is to compute the values of k
		
		std::size_t m = 0;
		for (std::size_t i = 0; i < 13; i++)
		{
			tpon = *at + alfa[i] * (*ah);
			
			for (std::size_t j = 0; j < x->size(); j++) xpon[j] = (*x)[j];
			
			for (std::size_t l = 0; l < i; l++)
			{
				++m;
                double beth = *ah*beta[m];

				for (std::size_t j = 0; j < x->size(); j++) xpon[j] += beth*k[l][j];
			}

			field(tpon, xpon, &dx);
			
			for (std::size_t j = 0; j < x->size(); j++) k[i][j] = dx[j];
		}
		
		if (*at < -0.76) { printf("%.10f, %.10f\n", dx[0], dx[1]); }
		
		//this is to compute the rk7 and rk8 predictions
		
		err = nor = 0.e0;
		
		for (std::size_t j = 0; j < x->size(); j++)
		{
			x7[j] = x8[j] = (*x)[j];
			
			for (std::size_t l = 0; l < 11; l++)
			{
                double kh = *ah*k[l][j];
				x7[j] += kh*c7[l];
				x8[j] += kh*c8[l];
			}

			x8[j] += *ah*(c8[11] * k[11][j] + c8[12] * k[12][j]);
			err += fabs(x8[j] - x7[j]);
			nor += fabs(x8[j]);
		}
		err /= x->size();
		
		//next lines compute the new time step h
		
		tol1 = tol * (1 + nor / 100);

		if (err < tol1) err = std::max(err, tol1 / 256);

		h1 = *ah;
		*ah *= 0.9 * pow(tol1 / err, 0.125);

		if (fabs(*ah) < hmin) *ah = hmin * SGN(*ah);

		if (fabs(*ah) > hmax) *ah = hmax * SGN(*ah);

	} while ((err >= tol1) && (fabs(*ah) > hmin));
	
	*at += h1;
	for (std::size_t j = 0; j<x->size(); j++) (*x)[j] = x8[j];

	return (err);
}
};

