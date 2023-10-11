#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
using namespace std;
ifstream fin;
ofstream fout;
#define PI 3.1415926535
int k = 100;
static long double f(const long double& x)
{
	return cos(x);
}
void e(double* x, double* y) {
	fout << "���������� ����������� ��� ����������:" << endl;
	for (int i = 0; i < k; i++) fout << "(" << x[i] << ";" << fabs(y[i] - f(x[i]))<< ")" << endl ;
	return;
}
void er(double* x, double* arr, int n) {
	fout << "������������� ����������� ��� ����������:" << endl;
	for (int i = 0; i < k; i++) {
		double e = 1;
		for (int j = 0; j < n; j++) {
			e *= (x[i] - arr[j]) / (j+1);
		}
		fout << "(" << x[i] << ";"<< fabs(e) << ")" << endl;
	}
	return;
}
void Interpol_Lagrange(double* x, double* arr, int n) {
	double* sum = new double[k];
	double an;
	for (int m = 0; m < k; m++) {
		sum[m] = 0;
		for (int i = 0; i < n; i++) {
			an = f(arr[i]);
			for (int j = 0; j < n; j++) {
				if (i != j) an *= (x[m] - arr[j]) / (arr[i] - arr[j]);
			}
			sum[m] += an;
		}
	}
	for (int i = 0; i < k; i++) {
		fout << "("<< x[i] << ";" << sum[i] << ")\n ";
	}
	e(x, sum);
	er(x, arr, n);
	return;

}
void  Interpol_Newton(double* x, double* arr, int n)
{
	double* sum = new double[k];
	for (int m = 0; m < k; m++) {
		sum[m] = 0;
		for (int j = 1; j <= n; j++) {
			double fj = 0;//
			for (int i = 0; i < j; i++) {
				double c = f(arr[i]);
				for (int l = 0; l < j; l++) if (l != i) c /= (arr[i] - arr[l]);
				fj += c;
			}
			double hsum = fj;
			for (int l = 0; l < j-1; l++) {
				hsum *= (x[m] - arr[l]);
			}
			sum[m] += hsum;
		}
	}
	for (int i = 0; i < k; i++) {
		fout << "(" << x[i] << ";" << sum[i] << ")\n ";
	}
	e(x, sum);
	er(x, arr, n);
	return;
}


int main() {
	double a, b; //������� ���������
	int n;// ��������� 
	srand(time(NULL));
	fin.open("input.txt");
	fout.open("output.txt");
	fin >> a >> b >> n;
	double* x = new double[k];
	double* arr1 = new double[n];
	double* arr2 = new double[n];
	double step = (b - a) / (n-1);
	arr2[0] = a; arr2[n - 1] = b;
	for (int i = 0; i < n; i++) {
		
		arr1[i] = a+step*i;// ������ ����������� �������� �����
		double d = (2 * i + 1) * PI / (2 * n + 2);
		if (i != 0 && i != n-1) {
			double h = 0.5 * ((b - a) * cos(d) + a + b);// ������ �� ������ ��������� ��������
			arr2[i] = h;
		}
	}
	//for (int i = 0; i < n; i++) {
	//	arr1[i] = a + step * i;// ������ ����������� �������� �����
	//}
	//for (int i = 1; i < n-1; i++) {

	//	
	//	double d = (2 * i + 1) * PI / (2 * n + 2);

	//		double h = 0.5 * ((b - a) * cos(d) + a + b);// ������ �� ������ ��������� ��������
	//		arr2[i] = h;
	//	
	//}
	for (int i = 0; i < k; i++) {
		x[i] = (double)rand() * (b - a) / RAND_MAX + a;
	}
	//for (int i = 0; i < n; i++) {
 //   	fout << "(" << arr1[i] << " ;" << f(arr1[i]) << ")" << endl;// ������ ����������� �������� �����
	//}
	for (int i = 0; i < n; i++) {
    	fout <<"(" <<arr2[i] << " ;" << f(arr2[i])<<")" << endl;// ������ ����������� �������� �����
	}
	/*fout << "���������������� ������� �������� � ����������� �����:\n";
	Interpol_Lagrange(x, arr1, n);*/
	fout << "\n\n ���������������� ������� �������� � ����� - ������ �������� ��������:\n";
	Interpol_Lagrange(x, arr2, n);
	/*fout << "���������������� ������� ������� � ����������� �����:\n";
	Interpol_Newton(x, arr1, n);
	fout << "\n\n ���������������� ������� ������� � ����� - ������ �������� ��������:\n";
	Interpol_Newton(x, arr2, n);*/
	fin.close();
	fout.close();
	return 0;
}