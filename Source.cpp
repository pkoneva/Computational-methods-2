#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
using namespace std;
ifstream fin;
ofstream fout;
#define PI 3.1415926535


/*----подинтегральная функция---*/
static double f(const double& x)
{
	return sin(sin(x))/*exp(-(x * x) / 2)*/;
}

static double f2(const double& x)
{
	return /*log(x)*/x*x;
}

/*----Интегрирование методом монте-карло---*/
static double Monte_Karlo(const double& left, const double& right, const double& M)
{
	double sum = 0;
	double x;
	for (int i = 0; i < M; i++) {
		x = (double)rand()* (right - left) / RAND_MAX + left;/*(double)((rand() / (right - left + 1)) + left);*/
		
		sum += f(x)*(right - left) / M;
	}

	const double Dsum = sum / (right - left);
	double D=0;
	for (int i = 0; i < M; i++) {
		x = (double)rand() * (right - left) / RAND_MAX + left;/*(double)((rand() / (right - left + 1)) + left);*/

		D += pow(f(x) - Dsum, 2)/(M-1);
	}
	D = sqrt(D);
	fout << "Дисперсия для метода Монте-Карло: " << D;
	fout << "\nОценка точности вычислений для метода Монте-Карло: " << 3 * /*D/M*/D/sqrt(M);
	return sum/**(right-left) / M*/;
	
}

/*----Интегрирование методом правых прямоугольников (составная формула) ---*/
static double Pr_pryamoyg(const double& left, const double& right, const double& h) {
	double sum = 0;
	double x = left;
	double len = (right - left) / h;
	for (int i = 1; i <= h; i++) {
		sum += f2(x + len * i) * len;
	}
	return sum;
}

int main() {
	double a, b;
	double a2, b2;
	double h, M;
	fin.open("input.txt");
	fout.open("output.txt");
	fin >> a >> b;
	fin >> a2 >> b2;
	fin >> h >> M;
	srand(time(NULL));

	fout <<"\nВычисление методом Монте-Карло: "<< Monte_Karlo(a, b, M);
	fout << "\nВычисление методом правых прямоугольников: " << Pr_pryamoyg(a2, b2, h);
	fin.close();
	fout.close();
	return 0;
}