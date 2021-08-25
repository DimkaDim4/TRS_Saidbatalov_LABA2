#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>

#define a0 0
#define a1 0
#define b0 1
#define b1 1

using namespace std;

double F(double x, double t)
{
	return (M_PI * t + 1) * cos(M_PI * x) - x / 2.;
}
double fi(double x)
{
	return sin(M_PI * x);
}
double psi_0(double t)
{
	return t;
}
double psi_1(double t)
{
	return -M_PI * exp(-M_PI * t);
}
double Solution(double x, double t)
{
	return exp(-M_PI * t) * sin(M_PI * x) + t * cos(M_PI * x) + t * x / 2;
}
double Error(double* y, double* y_sol, int N, double h, double tau)
{
	double maxdelta = -1.0;
	for (int j = 0; j <= N; j++)
	{
		for (int i = 0; i <= N; i++)
		{
			double delta = abs(y[j * (N + 1) + i] - y_sol[j * (N + 1) + i]);
			if (delta > maxdelta)
				maxdelta = delta;
		}
	}
	return maxdelta;
}
void Progonka(double* A, double* y, int j, int N)
{
	for (int i = 1; i < N; i++)
	{
		A[i + 1 * N] = A[i + 1 * N] - A[i - 1 + 0 * N] * (A[i - 1 + 2 * N] / A[i - 1 + 1 * N]);
		A[i + 3 * N] = A[i + 3 * N] - A[i - 1 + 3 * N] * (A[i - 1 + 2 * N] / A[i - 1 + 1 * N]);
		A[i - 1 + 2 * N] = 0;
	}

	y[N - 1 + N * j] = A[N - 1 + 3 * N] / A[N - 1 + 1 * N];

	for (int i = N - 2; i >= 0; i--)
		y[i + N * j] = (A[i + 3 * N] - A[i + 0 * N] * y[i + 1 + N * j]) / A[i + 1 * N];
}

//ЗАДАНИЕ 1
void KRS(int N, double h, double tau, double* y, double A)
{

	float start = clock();
	for (int i = 0; i <= N; i++)
		y[0 * (N + 1) + i] = fi(i * h);

	for (int n = 0; n < N; n++)
	{
		for (int i = 1; i < N; i++)

			y[(n + 1) * (N + 1) + i] = y[n * (N + 1) + i] + tau * F(i * h, n * tau) +
			A * (tau / (h * h)) * (y[n * (N + 1) + i + 1] - 2.0 * y[n * (N + 1) + i] + y[n * (N + 1) + i - 1]);

		y[(n + 1) * (N + 1) + 0] = (h * psi_0((n + 1) * tau) - b0 * y[(n + 1) * (N + 1) + 1]) / (a0 * h - b0);

		y[(n + 1) * (N + 1) + N] = (h * psi_1((n + 1) * tau) + b1 * y[(n + 1) * (N + 1) + N - 1]) / (a1 * h + b1);
	}
	cout << "\ntime: " << (clock() - start) / CLOCKS_PER_SEC << endl;
}

// ЗАДАНИЕ 2
void KRS_IMP(int N, double h, double tau, double* y, double a)
{

	for (int i = 0; i <= N; i++)
		y[0 * (N + 1) + i] = fi(i * h);
	float start = clock();

	for (int j = 0; j < N; j++)
	{
		double* A = new double[(N + 1) * 4];
		A[0 * (N + 1) + 0] = b0;
		A[1 * (N + 1) + 0] = a0 * h - b0;
		A[1 * (N + 1) + N] = a1 * h + b1;
		A[2 * (N + 1) + N - 1] = -b1;
		A[3 * (N + 1) + 0] = h * psi_0(tau * (j + 1));
		A[3 * (N + 1) + N] = h * psi_1(tau * (j + 1));
		for (int i = 1; i < N; i++)
		{
			A[0 * (N + 1) + i] = a * tau / (h * h);
			A[1 * (N + 1) + i] = -1 - (a * 2 * tau / (h * h));
			A[2 * (N + 1) + i - 1] = a * tau / ((h * h));
			A[3 * (N + 1) + i] = -y[j * (N + 1) + i] - tau * F(i * h, tau * j);
		}

		Progonka(A, y, j + 1, N + 1);
		delete[] A;
	}
	cout << "\ntime: " << (clock() - start) / CLOCKS_PER_SEC << endl;
}
void KRS_KN(int N, double h, double tau, double* y, double a)
{


	for (int i = 0; i <= N; i++)
		y[0 * (N + 1) + i] = fi(i * h);
	float start = clock();
	for (int j = 0; j < N; j++)
	{
		double* A = new double[(N + 1) * 4];
		A[0 * (N + 1) + 0] = b0;
		A[1 * (N + 1) + 0] = a0 * h - b0;
		A[1 * (N + 1) + N] = a1 * h + b1;
		A[2 * (N + 1) + N - 1] = -b1;
		A[3 * (N + 1) + 0] = h * psi_0(tau * (j + 1));
		A[3 * (N + 1) + N] = h * psi_1(tau * (j + 1));
		for (int i = 1; i < N; i++)
		{
			int k = j * (N + 1) + i;
			A[0 * (N + 1) + i] = a * tau / (2.0 * h * h);
			A[1 * (N + 1) + i] = -(1 + (a * tau / (h * h)));
			A[2 * (N + 1) + i - 1] = a * tau / (2.0 * (h * h));
			A[3 * (N + 1) + i] = -((1 - a * tau / (h * h)) * y[k] + a * tau / (2.0 * (h * h)) * (y[k - 1] + y[k + 1]) + tau * F(i * h, tau * j));
		}


		Progonka(A, y, j + 1, N + 1);

		delete[] A;
	}
	cout << "\ntime: " << (clock() - start) / CLOCKS_PER_SEC << endl;
	return;
}

// ЗАДАНИЕ 3
double K(double u)
{
	return cos(u);
}
void KRS_CONS(int N, double h, double tau, double* y, double a)
{

	double ke, kw;
	double* x = new double[N + 1];
	double* t = new double[N + 1];

	for (int i = 0; i <= N; i++)
	{
		y[0 * (N + 1) + i] = fi(i * h);
		x[i] = i * h;
		t[i] = i * tau;

	}

	float start = clock();

	for (int j = 0; j < N; j++)
	{

		double* A = new double[(N + 1) * 4];

		A[0 * (N + 1) + 0] = b0;
		A[1 * (N + 1) + 0] = a0 * h - b0;
		A[1 * (N + 1) + N] = a1 * h + b1;
		A[2 * (N + 1) + N - 1] = -b1;
		A[3 * (N + 1) + 0] = h * psi_0(tau * (j + 1));
		A[3 * (N + 1) + N] = h * psi_1(tau * (j + 1));



		for (int i = 1; i < N; i++)
		{
			int k = j * (N + 1) + i;

			kw = K((y[k - 1] + y[k]) / 2.);
			ke = K((y[k + 1] + y[k]) / 2.);

			A[0 * (N + 1) + i] = -a * kw / (h * h);
			A[1 * (N + 1) + i] = (1. / tau + (ke + kw) * (a / (h * h)));
			A[2 * (N + 1) + i - 1] = -a / ((h * h)) * ke;
			A[3 * (N + 1) + i] = ((1. / tau - (1. - a) * (ke + kw) / (h * h)) * y[k]
				+ (1. - a) * (y[k + 1] * ke - y[k - 1] * kw) / (h * h)
				+ cos(y[k]) * F(i * h, tau * j));
		}

		Progonka(A, y, j + 1, N + 1);
		delete[] A;


	}

	cout << "\ntime: " << (clock() - start) / CLOCKS_PER_SEC << endl;

}

int main()
{
	int N = 1000;
	double a = 0., b = 1.;
	double h = (b - a) / N;
	double tau = h * h / 2;
	double* y = new double[(N + 1) * (N + 1)];
	double* y_sol = new double[(N + 1) * (N + 1)];
	for (int j = 0; j <= N; j++)
	{
		for (int i = 0; i <= N; i++)
		{
			y_sol[j * (N + 1) + i] = Solution(i * h, j * tau);
		}
	}
	KRS(N, h, tau, y, 1 / M_PI);
	cout << "Error: " << Error(y, y_sol, N, h, tau) << endl;
	KRS_IMP(N, h, tau, y, 1 / M_PI);
	cout << "Error: " << Error(y, y_sol, N, h, tau) << endl;
	KRS_KN(N, h, tau, y, 1 / M_PI);
	cout << "Error: " << Error(y, y_sol, N, h, tau) << endl;
	KRS_CONS(N, h, tau, y, 1 / M_PI);
	cout << "Error: " << Error(y, y_sol, N, h, tau) << endl;

	return 0;
}
