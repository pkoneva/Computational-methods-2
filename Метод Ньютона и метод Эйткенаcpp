#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;
ifstream fin;
ofstream fout;
#define PI 3.1415926535
#define e1 0.0000001
#define e2 0.000000000001

static long double f(const long double& x)
{
	return cos(x)-x;
}

static long double f1(const long double& x)//1-я производная f
{
	return -sin(x)-1;
}

void M_Newton(long double& x) {
	long double xi = x - (f(x) / f1(x));
	int iter=1;
	while (fabs(x - xi) > e1) {
		x = xi;
		xi = x - f(x) / f1(x);
		iter++;
	}
	fout << "Приближение №1 методом Ньютона: " << setprecision(6) << xi << "\nКоличество итераций: "<<iter;
	while (fabs(x - xi) > e2) {
		x = xi;
		xi = x - f(x) / f1(x);
		iter++;
	}
	fout << "\nПриближение №2 методом Ньютона: " << setprecision(12) <<xi << "\nКоличество итераций: " << iter;
	return;
}
void M_Aitken(long double& x0) {
	long double x1, x2, dx1, dx2;
	long double x;
	int iter = 0;
	do {
		iter++;
		x1 = cos(x0);
		x2 = cos(x1);
		dx1 = x1 - x0;
		dx2 = x2 - x1;
		x = x2 - dx2 * dx2 / dx1;
		x0 = x1;
		

	} while (fabs(x - x2) > e1);
	fout << "\nПриближение №1 методом Эйткена: " << setprecision(6) << x << "\nКоличество итераций: " << iter;
	while (fabs(x - x2) > e2) {
		iter++;
		x1 = cos(x0);
		x2 = cos(x1);
		dx1 = x2 - x0;
		dx2 = x2 - x1;
		x = x2 - pow((dx2 / dx1), 2);
		x0 = x1;
	}
	fout << "\nПриближение №2 методом Эйткена: " << setprecision(12) << x << "\nКоличество итераций: " << iter;
	return;
}


int main() {

	fin.open("input.txt");
	fout.open("output.txt");
	srand(time(NULL));
	long double x = (long double)rand() / RAND_MAX; //нормировка, чтобы значение попало в интервал 0,1
	M_Newton(x);
	fout << "\n\n";
	x = (long double)rand() / RAND_MAX;
	M_Aitken(x);
	fin.close();
	fout.close();
	return 0;
}
