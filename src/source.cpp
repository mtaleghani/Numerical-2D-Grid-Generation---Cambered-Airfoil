/* CFD II Project 2 Code
   Version Number: 4
   By Seyed MohammadAmin Taleghani
   Sharif University of Technology - Aerospace Engineering Department */

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
using namespace std;

#define imax 101 /* Number of the Nodes in Xsi-dir (SHOULD BE ALWAYS ODD)*/
#define jmax 81 /* Number of the Nodes in Etha-dir */
double err = 8e-6; /* defines the convergence critetion for iterations */

							// IMAX + 1 points for periodic BC
double x_jac[imax][jmax], x_PSOR[imax + 1][jmax], x_XiSweep[imax][jmax], x_EthaSweep[imax][jmax], x_adi[imax][jmax];
double x_sor[imax][jmax], x_new[imax + 1][jmax], T_anal[imax][jmax], diffx[imax][jmax] = { 0.0 }, diffy[imax][jmax] = { 0.0 };
// IMAX + 1 points for periodic BC
double y_jac[imax][jmax], y_PSOR[imax + 1][jmax], y_XiSweep[imax][jmax], y_EthaSweep[imax][jmax], y_adi[imax][jmax];
double y_sor[imax][jmax], y_new[imax + 1][jmax];

int iter_jac = 0, iter_gauss = 0, iter_rsweep = 0, iter_thsweep = 0, iter_adi_rth = 0, iter_adi_thr = 0, iter_sor = 0;

/* Function Declarations */

void TecplotOutput(const char* Name, double x[imax][jmax], double y[imax][jmax])
{
	ofstream MyFile(Name);

	/* Tecplot Output*/

	/* Header */
	MyFile << "Variables = x , y" << "\n";

	/* Header */

	/* Data */
	MyFile << "Zone T = \"Test\" , I = " << imax << " , J = " << jmax << "\n";
	for (int j = 0; j < jmax; j++)
	{
		for (int i = 0; i < imax; i++)
		{
			MyFile << x[i][j];
			MyFile << "		";
			MyFile << y[i][j];
			MyFile << "\n";
		}
	}
	/* Data */

	MyFile.close();
}

double NewtonRaphsDelta(double B)
{
	double err = 1e-5;
	double x = 0.5;

	do
	{
		x -= (sinh(x) / x - B) / (cosh(x) / x - sinh(x) / (x * x));		// x = x - (y(x)/y'(x))
	} while (sinh(x) / x - B >= err || sinh(x) / x - B <= -err);

	return x;
}

void Vink(double* x, double dsa, double dsb, double smax, int imaxu)
{
	double A = sqrt(dsb / dsa);
	double B = smax / ((imaxu - 1.0) * sqrt(dsa * dsb));
	double delta = NewtonRaphsDelta(B), eta;
	double* u;
	u = new double[imaxu];

	for (int i = 0; i < imaxu; i++)
	{
		eta = (i) / (imaxu - 1.0);
		u[i] = 0.5 + tanh(delta * (eta - 0.5)) / (2 * tanh(delta / 2));
		x[i] = smax * u[i] / (A + (1 - A) * u[i]);
	}
}

double maxvec(double arr[], int n)
{
	int i;

	// Initialize maximum element
	double max = arr[0];

	// Traverse array elements
	// from second and compare
	// every element with current max
	for (i = 1; i < n; i++)
		if (arr[i] > max)
			max = arr[i];

	return max;
}

double maxarr(double arr[imax][jmax], int n, int m)
{
	int i, j;

	// Initialize maximum element
	double max = arr[0][0];

	// Traverse array elements
	// from second and compare
	// every element with current max
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			if (arr[i][j] > max)
			{
				max = arr[i][j];
			}
		}
	}
	return max;
}

double* thomas(double A[][imax - 2], double B[])
{
	double ratio; static double X[imax - 2];

	for (int i = 1; i < imax - 2; i++)
	{
		ratio = A[i][i - 1] / A[i - 1][i - 1];
		A[i][i] = A[i][i] - ratio * A[i - 1][i];
		B[i] = B[i] - ratio * B[i - 1];
	}

	X[imax - 2 - 1] = B[imax - 2 - 1] / A[imax - 2 - 1][imax - 2 - 1];

	for (int i = imax - 2 - 2; i >= 0; i--)
	{
		X[i] = (B[i] - A[i][i + 1] * X[i + 1]) / A[i][i];
	}

	return X;
}

void NacaSymm(double* x_vink, double x[imax][jmax], double y[imax][jmax], int iminu, int imaxu, double t)
{
	for (int i = iminu; i < imaxu; i++)
	{
		/*x[i][0] = -0.5 + i / (imaxu - 1.0);*/  //this was for the uniform grid on the airfoil
		x[i][0] = -0.5 + x_vink[i];
		double x = x_vink[i];
		y[i][0] = 5.0 * t * (0.2969 * sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x * x * x);
	}

	y[imaxu - 1][0] = 0.0;

	for (int i = imaxu; i < imax; i++)
	{
		x[i][0] = x[2 * imaxu - 2 - i][0];
		y[i][0] = -y[2 * imaxu - 2 - i][0];
	}

	x[imax][0] = x[1][0];
	y[imax][0] = y[1][0];
}

void NacaCamb(double* x_vink, double x[imax][jmax], double y[imax][jmax], int iminu, int imaxu, double m, double p, double t)
{
	const int imu = 0.5 * (imax - 1.0) + 1;

	double yt[imu], xc[imu], yc[imu], dyc[imu], th[imu];

	for (int i = iminu; i < imaxu; i++)
	{

		xc[i] = x_vink[i];
		double x = x_vink[i];
		yt[i] = 5.0 * t * (0.2969 * sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x * x * x);
		if (x <= p)
		{
			yc[i] = m / (p * p) * (2. * p * x - x * x);
			dyc[i] = 2. * m / (p * p) * (p - x);
			th[i] = atan(dyc[i]);
		}

		if (p < x)
		{
			yc[i] = m / pow(1. - p, 2) * ((1. - 2. * p) + 2. * p * x - x * x);
			dyc[i] = 2. * m / pow(1. - p, 2) * (p - x);
			th[i] = atan(dyc[i]);
		}

	}


	for (int i = 1; i < imaxu; i++)
	{
		x[i][0] = xc[i] - yt[i] * sin(th[i]) - 0.5;
		y[i][0] = yc[i] + yt[i] * cos(th[i]);

		x[imax - 1 - i][0] = xc[i] + yt[i] * sin(th[i]) - 0.5;
		y[imax - 1 - i][0] = yc[i] - yt[i] * cos(th[i]);
	}


	x[0][0] = -0.5;
	y[0][0] = 0.0;

	x[imaxu - 1][0] = 0.5;
	y[imaxu - 1][0] = 0.0;

	x[imax - 1][0] = -0.5;
	y[imax - 1][0] = 0.0;

	x[imax][0] = x[1][0];
	y[imax][0] = y[1][0];

}

void Cylind(double* x_vink, double x[imax][jmax], double y[imax][jmax], double r)
{

	for (int i = 0; i <= (imax - 1) / 2; i++)
	{
		x[i][0] = -r + x_vink[i];
		y[i][0] = sqrt(r * r - x[i][0] * x[i][0]);
	}

	for (int i = 1; i < (imax - 1) / 2; i++)
	{
		x[(imax - 1) / 2 + i][0] = x[(imax - 1) / 2 - i][0];
		y[(imax - 1) / 2 + i][0] = -y[(imax - 1) / 2 - i][0];
	}

	x[imax - 1][0] = x[0][0];
	y[imax - 1][0] = y[0][0];

	x[imax][0] = x[1][0];
	y[imax][0] = y[1][0];
}

void present(double* x_vink, double x[imax][jmax], double y[imax][jmax], double r)
{
	for (int i = 0; i <= (imax - 1) / 4; i++)
	{
		x[i][0] = -r * cos(x_vink[i]);
		y[i][0] = r + r * sin(x_vink[i]);
	}

	for (int i = 0; i <= (imax - 1) / 4; i++)
	{
		x[(imax - 1) / 4 + i][0] = r - r * sin(-x_vink[i]);
		y[(imax - 1) / 4 + i][0] = r * cos(-x_vink[i]);
	}

	for (int i = 0; i <= (imax - 1) / 4; i++)
	{
		x[(imax - 1) / 2 + i][0] = r * cos(-x_vink[i]);
		y[(imax - 1) / 2 + i][0] = -r + r * sin(-x_vink[i]);
	}

	for (int i = 0; i <= (imax - 1) / 4; i++)
	{
		x[3 * (imax - 1) / 4 + i][0] = -r - r * sin(x_vink[i]);
		y[3 * (imax - 1) / 4 + i][0] = -r * cos(x_vink[i]);
	}


	x[imax - 1][0] = x[0][0];
	y[imax - 1][0] = y[0][0];

	x[imax][0] = x[1][0];
	y[imax][0] = y[1][0];

}


void Initialize(double x[imax][jmax], double y[imax][jmax], bool algebric)
{
	// Discretize Outer Circle
	double r = 12.0;

	for (int i = 0; i < imax + 1; i++)
	{
		x[i][jmax - 1] = r * cos(M_PI * (1.0 - 2.0 * i / (imax - 1.0)));
		y[i][jmax - 1] = r * sin(M_PI * (1.0 - 2.0 * i / (imax - 1.0)));
	}

	if (algebric)
	{
		// Interpolate the mesh inside the domain.
		double dx, dy;

		for (int j = 1; j < jmax - 1; j++)
		{
			for (int i = 0; i < imax; i++)
			{
				dx = (x[i][jmax - 1] - x[i][0]) / (jmax - 1.0);
				dy = (y[i][jmax - 1] - y[i][0]) / (jmax - 1.0);
				x[i][j] = x[i][0] + j * dx;
				y[i][j] = y[i][0] + j * dy;
			}

			x[imax][j] = x[1][j];
			y[imax][j] = y[1][j];
		}
	}

	if (!algebric)
	{
		double Guessx = 0.0, Guessy = 0.0;

		for (int j = 1; j < jmax - 1; j++)
		{
			for (int i = 0; i < imax; i++)
			{
				x[i][j] = Guessx;
				y[i][j] = Guessy;
			}

			x[imax][j] = x[1][0];
			y[imax][j] = y[1][0];
		}
	}
}

void SolvePSOR(double x[][jmax], double y[][jmax], double ds, bool Auto, int iter_start, double Q1, double P1, double a, double b, double omega, double omegaPQ)
{
	double xsi0 = 0.0, eta0 = 0.0, dxsi = 1.0, deta = 1.0;   /*dxsi = 1.0 / (imax - 1), deta = 1.0 / (jmax - 1);*/
	//double xsi[imax], eta[jmax];
	double L_infty = 0, L_inftx = 0, sum = 0.0;

	int i, j;

	double Alpha, Beta, Gamma, xxsi, XXSI, yxsi, YXSI, xeta, XETA, yeta, YETA, xxsi2, xxsieta, xeta2, yxsi2, yxsieta, yeta2, sourcex[imax + 1] = { 0.0 }, sourcey[imax + 1] = { 0.0 };
	double J, R1, R2, P1vec[imax + 1] = { 0.0 }, Q1vec[imax + 1] = { 0.0 };

	// Tic
	std::clock_t c_start = std::clock();

	int iter = 0;

	for (i = 0; i < imax + 1; i++)
	{
		P1vec[i] = P1;
		Q1vec[i] = Q1;
	}

	while (1)
	{
		//calculate P1 , Q1 for Boundary j=0

		if (Auto)
		{
			j = 0;

			for (i = 1; i < imax; i++)
			{
				xxsi = (x[i + 1][j] - x[i - 1][j]) / (2.0 * dxsi);
				yxsi = (y[i + 1][j] - y[i - 1][j]) / (2.0 * dxsi);

				xeta = -(ds / deta) * yxsi / sqrt(pow(xxsi, 2) + pow(yxsi, 2)); // signs as per jozve
				yeta = (ds / deta) * xxsi / sqrt(pow(xxsi, 2) + pow(yxsi, 2));

				J = xxsi * yeta - xeta * yxsi;
				Alpha = (xeta * xeta + yeta * yeta);
				Beta = (xxsi * xeta + yxsi * yeta);
				Gamma = (xxsi * xxsi + yxsi * yxsi);

				xeta2 = (-7. * x[i][j] + 8. * x[i][j + 1] - x[i][j + 2]) / (2. * deta * deta) - 3. * xeta / deta;
				yeta2 = (-7. * y[i][j] + 8. * y[i][j + 1] - y[i][j + 2]) / (2. * deta * deta) - 3. * yeta / deta;

				xxsieta = (x[i + 1][j + 1] - x[i + 1][j] + x[i - 1][j] - x[i - 1][j + 1]) / (2. * dxsi * deta); // O( dxsi^2 , deta )
				yxsieta = (y[i + 1][j + 1] - y[i + 1][j] + y[i - 1][j] - y[i - 1][j + 1]) / (2. * dxsi * deta); // O( dxsi^2 , deta )

				xxsi2 = (x[i + 1][j] - 2. * x[i][j] + x[i - 1][j]) / (dxsi * dxsi);
				yxsi2 = (y[i + 1][j] - 2. * y[i][j] + y[i - 1][j]) / (dxsi * dxsi);

				R1 = -(Alpha * xxsi2 - 2. * Beta * xxsieta + Gamma * xeta2) / (J * J);
				R2 = -(Alpha * yxsi2 - 2. * Beta * yxsieta + Gamma * yeta2) / (J * J);

				P1vec[i] = omegaPQ * -(+yeta * R1 - xeta * R2) / J + (1. - omegaPQ) * P1vec[i]; // As per Steger & Sorenson
				Q1vec[i] = omegaPQ * -(-yxsi * R1 + xxsi * R2) / J + (1. - omegaPQ) * Q1vec[i]; //
			}

			if (i == imax - 1)
			{
				P1vec[0] = P1vec[imax - 1];
				Q1vec[0] = Q1vec[imax - 1];
			}

			if (i == 1)
			{
				P1vec[imax] = P1vec[1];
				Q1vec[imax] = Q1vec[1];
			}
		}

		for (i = 1; i < imax; i++)
		{
			for (j = 1; j < jmax - 1; j++)
			{
				//calculate derivatives J, alpha , beta , gamma

				xeta = (x[i][j + 1] - x[i][j - 1]) / (2.0 * deta);
				yeta = (y[i][j + 1] - y[i][j - 1]) / (2.0 * deta);
				xxsi = (x[i + 1][j] - x[i - 1][j]) / (2.0 * dxsi);
				yxsi = (y[i + 1][j] - y[i - 1][j]) / (2.0 * dxsi);

				J = xxsi * yeta - xeta * yxsi;
				Alpha = (xeta * xeta + yeta * yeta);
				Beta = (xxsi * xeta + yxsi * yeta);
				Gamma = (xxsi * xxsi + yxsi * yxsi);

				//calculcate Q , P control function

				if (P1vec[i] == 0)
				{
					XXSI = (x[i + 1][j] - x[i - 1][j]) / (2.0 * dxsi);
					YXSI = (y[i + 1][j] - y[i - 1][j]) / (2.0 * dxsi);
				}
				if (P1vec[i] > 0)   // >=
				{
					XXSI = (x[i + 1][j] - x[i][j]) / (dxsi);
					YXSI = (y[i + 1][j] - y[i][j]) / (dxsi);
				}
				if (P1vec[i] < 0)
				{
					XXSI = (x[i][j] - x[i - 1][j]) / (dxsi);
					YXSI = (y[i][j] - y[i - 1][j]) / (dxsi);
				}

				if (Q1vec[i] == 0)
				{
					XETA = (x[i][j + 1] - x[i][j - 1]) / (2.0 * deta);
					YETA = (y[i][j + 1] - y[i][j - 1]) / (2.0 * deta);
				}
				if (Q1vec[i] > 0) // >=
				{
					XETA = (x[i][j + 1] - x[i][j]) / (deta);
					YETA = (y[i][j + 1] - y[i][j]) / (deta);
				}
				if (Q1vec[i] < 0)
				{
					XETA = (x[i][j] - x[i][j - 1]) / (deta);
					YETA = (y[i][j] - y[i][j - 1]) / (deta);
				}

				//calculate new x,y
				x_new[i][j] = ((x[i + 1][j] + x[i - 1][j]) * Alpha + (x[i][j + 1] + x[i][j - 1]) * Gamma - (x[i + 1][j + 1] - x[i + 1][j - 1] + x[i - 1][j - 1] - x[i - 1][j + 1]) * Beta * 0.5 + J * J * (P1vec[i] * exp(-a * deta * j) * XXSI + Q1vec[i] * exp(-b * deta * j) * XETA)) / (2. * (Alpha + Gamma));
				y_new[i][j] = ((y[i + 1][j] + y[i - 1][j]) * Alpha + (y[i][j + 1] + y[i][j - 1]) * Gamma - (y[i + 1][j + 1] - y[i + 1][j - 1] + y[i - 1][j - 1] - y[i - 1][j + 1]) * Beta * 0.5 + J * J * (P1vec[i] * exp(-a * deta * j) * YXSI + Q1vec[i] * exp(-b * deta * j) * YETA)) / (2. * (Alpha + Gamma));

				x_new[i][j] = omega * x_new[i][j] + (1. - omega) * x[i][j];
				y_new[i][j] = omega * y_new[i][j] + (1. - omega) * y[i][j];

				diffx[i][j] = abs((x_new[i][j] - x[i][j]));
				diffy[i][j] = abs((y_new[i][j] - y[i][j]));

				x[i][j] = x_new[i][j];
				y[i][j] = y_new[i][j];

				if (i == imax - 1)
				{
					x[0][j] = x[imax - 1][j];
					y[0][j] = y[imax - 1][j];
				}

				if (i == 1)
				{
					x[imax][j] = x[1][j];
					y[imax][j] = y[1][j];
				}
			}
		}

		cout << "Iter is : " << iter << "\n\n";
		L_inftx = maxarr(diffx, imax + 1, jmax + 1);
		cout << "L_inftx is : " << L_inftx << "\n\n";
		L_infty = maxarr(diffy, imax + 1, jmax + 1);
		cout << "L_infty is : " << L_infty << "\n\n";

		if (L_infty < err && L_inftx < err)
		{
			cout << "Iter Guass is : " << iter;

			break;
		}

		iter++;
	}

	// Toc
	//std::clock_t c_end = std::clock();
	//long double cpu_time_gauss = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
}

void SolveLSOR(double x[][jmax], double y[][jmax], double ds, bool Auto, int iter_start, double Q1, double P1, double a, double b, double omega, double omegaPQ)
{
	double xsi0 = 0.0, eta0 = 0.0, dxsi = 1.0, deta = 1.0;  // dxsi = 1.0 / (imax - 1), deta = 1.0 / (jmax - 1);
	//double xsi[imax], eta[jmax];
	double L_infty = 0, L_inftx = 0, sum = 0.0;

	int i, j;

	double Alpha, Beta, Gamma, xxsi, XXSI, yxsi, YXSI, xeta, XETA, yeta, YETA, xxsi2, xxsieta, xeta2, yxsi2, yxsieta, yeta2, sourcex[imax + 1] = { 0.0 }, sourcey[imax + 1] = { 0.0 };
	double J, R1, R2, P1vec[imax + 1] = { 0.0 }, Q1vec[imax + 1] = { 0.0 };

	//////// Defines A, B and X matrices for A*X=B;
	double A[imax - 2][imax - 2] = { 0.0 }, B[imax - 2];
	double* X;
	double* Y;

	double x_old[imax][jmax] = { 0.0 },
		y_old[imax][jmax] = { 0.0 };

	// Tic
	std::clock_t c_start = std::clock();

	int iter = 0;

	for (i = 0; i < imax + 1; i++)
	{
		P1vec[i] = P1;
		Q1vec[i] = Q1;
	}

	while (1)
	{
		//calculate P1 , Q1 for Boundary j=0

		if (Auto)
		{
			j = 0;

			for (i = 1; i < imax; i++)
			{
				xxsi = (x[i + 1][j] - x[i - 1][j]) / (2.0 * dxsi);
				yxsi = (y[i + 1][j] - y[i - 1][j]) / (2.0 * dxsi);

				xeta = -(ds / deta) * yxsi / sqrt(pow(xxsi, 2) + pow(yxsi, 2)); // signs as per jozve
				yeta = (ds / deta) * xxsi / sqrt(pow(xxsi, 2) + pow(yxsi, 2));

				J = xxsi * yeta - xeta * yxsi;
				Alpha = (xeta * xeta + yeta * yeta);
				Beta = (xxsi * xeta + yxsi * yeta);
				Gamma = (xxsi * xxsi + yxsi * yxsi);

				xeta2 = (-7. * x[i][j] + 8. * x[i][j + 1] - x[i][j + 2]) / (2. * deta * deta) - 3. * xeta / deta;
				yeta2 = (-7. * y[i][j] + 8. * y[i][j + 1] - y[i][j + 2]) / (2. * deta * deta) - 3. * yeta / deta;

				xxsieta = (x[i + 1][j + 1] - x[i + 1][j] + x[i - 1][j] - x[i - 1][j + 1]) / (2. * dxsi * deta); // O( dxsi^2 , deta )
				yxsieta = (y[i + 1][j + 1] - y[i + 1][j] + y[i - 1][j] - y[i - 1][j + 1]) / (2. * dxsi * deta); // O( dxsi^2 , deta )

				xxsi2 = (x[i + 1][j] - 2. * x[i][j] + x[i - 1][j]) / (dxsi * dxsi);
				yxsi2 = (y[i + 1][j] - 2. * y[i][j] + y[i - 1][j]) / (dxsi * dxsi);

				R1 = -(Alpha * xxsi2 - 2. * Beta * xxsieta + Gamma * xeta2) / (J * J);
				R2 = -(Alpha * yxsi2 - 2. * Beta * yxsieta + Gamma * yeta2) / (J * J);

				P1vec[i] = omegaPQ * -(+yeta * R1 - xeta * R2) / J + (1. - omegaPQ) * P1vec[i]; // As per Steger & Sorenson
				Q1vec[i] = omegaPQ * -(-yxsi * R1 + xxsi * R2) / J + (1. - omegaPQ) * Q1vec[i]; //
			}

			if (i == imax - 1)
			{
				P1vec[0] = P1vec[imax - 1];
				Q1vec[0] = Q1vec[imax - 1];
			}

			if (i == 1)
			{
				P1vec[imax] = P1vec[1];
				Q1vec[imax] = Q1vec[1];
			}
		}
		// end of if

		////////////// Sets T_old = T_rsweep

		for (i = 0; i < imax; i++)
		{
			for (j = 0; j < jmax; j++)
			{
				x_old[i][j] = x[i][j];
				y_old[i][j] = y[i][j];
			}
		}

		for (j = 1; j < jmax - 1; j++)
		{
			//calculate derivatives J, alpha , beta , gamma

			xeta = (x[i][j + 1] - x[i][j - 1]) / (2.0 * deta);
			yeta = (y[i][j + 1] - y[i][j - 1]) / (2.0 * deta);
			xxsi = (x[i + 1][j] - x[i - 1][j]) / (2.0 * dxsi);
			yxsi = (y[i + 1][j] - y[i - 1][j]) / (2.0 * dxsi);

			J = xxsi * yeta - xeta * yxsi;
			Alpha = (xeta * xeta + yeta * yeta);
			Beta = (xxsi * xeta + yxsi * yeta);
			Gamma = (xxsi * xxsi + yxsi * yxsi);

			//calculcate Q , P control function

			if (P1vec[i] == 0)
			{
				XXSI = (x[i + 1][j] - x[i - 1][j]) / (2.0 * dxsi);
				YXSI = (y[i + 1][j] - y[i - 1][j]) / (2.0 * dxsi);
			}
			if (P1vec[i] > 0)   // >=
			{
				XXSI = (x[i + 1][j] - x[i][j]) / (dxsi);
				YXSI = (y[i + 1][j] - y[i][j]) / (dxsi);
			}
			if (P1vec[i] < 0)
			{
				XXSI = (x[i][j] - x[i - 1][j]) / (dxsi);
				YXSI = (y[i][j] - y[i - 1][j]) / (dxsi);
			}

			if (Q1vec[i] == 0)
			{
				XETA = (x[i][j + 1] - x[i][j - 1]) / (2.0 * deta);
				YETA = (y[i][j + 1] - y[i][j - 1]) / (2.0 * deta);
			}
			if (Q1vec[i] > 0) // >=
			{
				XETA = (x[i][j + 1] - x[i][j]) / (deta);
				YETA = (y[i][j + 1] - y[i][j]) / (deta);
			}
			if (Q1vec[i] < 0)
			{
				XETA = (x[i][j] - x[i][j - 1]) / (deta);
				YETA = (y[i][j] - y[i][j - 1]) / (deta);
			}

			//////////// Fills A, B and X matrices for r-sweep (the 1st and last row formulas are different and are filled outside the loop);

			A[0][0] = -2.0 * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta));
			A[0][1] = omega * Alpha / (dxsi * dxsi);
			B[0] = -2. * (1.0 - omega) * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta)) * x_old[i][j] + omega * (x_old[i + 1][j + 1] - x[i + 1][j - 1] + x[i - 1][j - 1] - x_old[i - 1][j + 1]) * Beta * 0.5 / (dxsi * deta) - omega * (x_old[i][j + 1] + x[i][j - 1]) * Gamma / (deta * deta) - Alpha / (dxsi * dxsi) * x[imax - 1][j];

			for (i = 1; i < imax - 2 - 1; i++)
			{
				A[i][i - 1] = omega * Alpha / (dxsi * dxsi);
				A[i][i] = -2.0 * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta));
				A[i][i + 1] = omega * Alpha / (dxsi * dxsi);
				B[i] = -2. * (1.0 - omega) * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta)) * x_old[i][j] + omega * (x_old[i + 1][j + 1] - x[i + 1][j - 1] + x[i - 1][j - 1] - x_old[i - 1][j + 1]) * Beta * 0.5 / (dxsi * deta) - omega * (x_old[i][j + 1] + x[i][j - 1]) * Gamma / (deta * deta);
			}

			A[imax - 2 - 1][imax - 2 - 2] = omega * Alpha / (dxsi * dxsi);
			A[imax - 2 - 1][imax - 2 - 1] = -2.0 * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta));
			B[imax - 2 - 1] = -2. * (1.0 - omega) * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta)) * x_old[i][j] + omega * (x_old[i + 1][j + 1] - x[i + 1][j - 1] + x[i - 1][j - 1] - x_old[i - 1][j + 1]) * Beta * 0.5 / (dxsi * deta) - omega * (x_old[i][j + 1] + x[i][j - 1]) * Gamma / (deta * deta) - Alpha / (dxsi * dxsi) * x[0][j];

			X = thomas(A, B);

			for (i = 0; i < imax - 2; i++)
			{
				x[i + 1][j] = *(X + i);
			}

			x[imax - 1][j] = x[0][j];

			///////////// Fills A, B and X matrices for r-sweep (the 1st and last row formulas are different and are filled outside the loop);

			A[0][0] = -2.0 * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta));
			A[0][1] = omega * Alpha / (dxsi * dxsi);
			B[0] = -2. * (1.0 - omega) * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta)) * y_old[i][j] + omega * (y_old[i + 1][j + 1] - y[i + 1][j - 1] + y[i - 1][j - 1] - y_old[i - 1][j + 1]) * Beta * 0.5 / (dxsi * deta) - omega * (y_old[i][j + 1] + y[i][j - 1]) * Gamma / (deta * deta) - Alpha / (dxsi * dxsi) * y[imax - 1][j];

			for (i = 1; i < imax - 2 - 1; i++)
			{
				A[i][i - 1] = omega * Alpha / (dxsi * dxsi);
				A[i][i] = -2.0 * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta));
				A[i][i + 1] = omega * Alpha / (dxsi * dxsi);
				B[i] = -2. * (1.0 - omega) * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta)) * y_old[i][j] + omega * (y_old[i + 1][j + 1] - y[i + 1][j - 1] + y[i - 1][j - 1] - y_old[i - 1][j + 1]) * Beta * 0.5 / (dxsi * deta) - omega * (y_old[i][j + 1] + y[i][j - 1]) * Gamma / (deta * deta);
			}

			A[imax - 2 - 1][imax - 2 - 2] = omega * Alpha / (dxsi * dxsi);
			A[imax - 2 - 1][imax - 2 - 1] = -2.0 * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta));
			B[imax - 2 - 1] = -2. * (1.0 - omega) * (Alpha / (dxsi * dxsi) + Gamma / (deta * deta)) * y_old[i][j] + omega * (y_old[i + 1][j + 1] - y[i + 1][j - 1] + y[i - 1][j - 1] - y_old[i - 1][j + 1]) * Beta * 0.5 / (dxsi * deta) - omega * (y_old[i][j + 1] + y[i][j - 1]) * Gamma / (deta * deta) - Alpha / (dxsi * dxsi) * y[0][j];

			Y = thomas(A, B);

			for (i = 0; i < imax - 2; i++)
			{
				y[i + 1][j] = *(Y + i);
			}

			y[imax - 1][j] = y[0][j];

			for (i = 0; i < imax; i++)
			{
				for (j = 0; j < jmax; j++)
				{
					diffx[i][j] = abs((x[i][j] - x_old[i][j]));
					diffy[i][j] = abs((y[i][j] - y_old[i][j]));
				}
			}

			cout << "Iter is : " << iter << "\n\n";
			L_inftx = maxarr(diffx, imax + 1, jmax + 1);
			cout << "L_inftx is : " << L_inftx << "\n\n";
			L_infty = maxarr(diffy, imax + 1, jmax + 1);
			cout << "L_infty is : " << L_infty << "\n\n";

			if (L_infty < err && L_inftx < err)
			{
				cout << "Iter Guass is : " << iter;

				break;
			}
		}

		iter++;
	}
}

int main()
{

	// Symm & Cambered Airoils
	// Determine the index of trailing edge.
	int iminu = 0;
	int imaxu = 0.5 * (imax - 1) + 1;

	// Vinokur's streching function.
	double* x_vink;
	x_vink = new double[imaxu];


	//// Symm Airfoil
	//// Airfoil thickness wrt unit chord.
	//double t = 0.12;
	//double dsa = 0.01, dsb = 0.02, smax = 1.0;
	//Vink(x_vink, dsa, dsb, smax, imaxu);

	//// Discretize Airfoil's Surface.
	//NacaSymm(x_vink, x_PSOR, y_PSOR, iminu, imaxu, t);
	////


	//// Cambered Airfoil
	//// maximum camber (100 m is the first of the four digits)
	//double m = 4. / 100;
	//// location of maximum camber (10 p is the second digit in the NACA xxxx description)
	//double p = 4. / 10;
	//// thickness wrt unit chord
	//double t = .12;
	//double dsa = 0.01, dsb = 0.01, smax = 1.0;
	//Vink(x_vink, dsa, dsb, smax, imaxu);
	//NacaCamb(x_vink, x_PSOR, y_PSOR, iminu, imaxu, m, p, t);


	// Cylinder
	// Inner Cylinder's Radius.
	double r = 1.;
	double dsa = 0.0005, dsb = 0.0005, smax = 2. * r;
	Vink(x_vink, dsa, dsb, smax, (imax - 1) / 2 + 1);

	// Discretize Cylinder's Surface
	Cylind(x_vink, x_PSOR, y_PSOR, r);
	//


	//// Cylinder
	//// Inner Cylinder's Radius.
	//double r = 1.;
	//double dsa = M_PI / ((imax - 1) / 4.), dsb = M_PI / ((imax - 1) / 4.), smax = M_PI;
	//Vink(x_vink, dsa, dsb, smax, (imax - 1) / 4 + 1);

	////Discretize Inner Surface
	//present(x_vink, x_PSOR, y_PSOR, r);




	// Algebric grid inside the domain.
	bool algebric = true;
	Initialize(x_PSOR, y_PSOR, algebric);

	// Print initialized grid.
	TecplotOutput("Initialized.plt", x_PSOR, y_PSOR);

	double dsxi = 1.0 / (imax - 1), deta = 1.0 / (jmax - 1);

	// Solver (Gauss-Sidel)
	double ds = -0.005;
	double P1 = -0.00000000000,  // +- 0.1
		Q1 = -0.00000000000,  //  -0.00005
		a = 1.0, b = 1.0, omegaPQ = 1e-7, omega = 1.5;
	bool Auto = true;
	int iter_start = 0;

	SolvePSOR(x_PSOR, y_PSOR, ds, Auto, iter_start, Q1, P1, a, b, omega, omegaPQ);
	//SolveLSOR(x_PSOR, y_PSOR, ds, Auto, iter_start, Q1, P1, a, b, omega, omegaPQ);

	TecplotOutput("Final.plt", x_PSOR, y_PSOR);

	printf("\n\nPress Enter to terminate the program.\n");
	cin.get();

	return 0;
}
