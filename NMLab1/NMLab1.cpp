// NM1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <math.h>
#include <cstdio>
#include <fstream>
#include <iostream>

typedef float form;
typedef float form2;

class Matrix
{
public:
	Matrix *M = this;
	form **al;
	form **au;
	form *di;

	int n;
	int m;
	int p;

	void Decompose()
	{
		for (int i = 0; i < n; i++)
		{
			form2 sumL;
			form2 sumU;
			for (int j = i-p, jl=0; j<=i; j++, jl++)
			{
				if (j < 0) continue;
				sumL = 0;
				sumU = 0;

				for (int k = 0,ku=i-j; k < jl; k++, ku++)
				{
					
					//sumL = L(i, k)*U(k, j);
					sumL += al[i][k] * au[j][ku];
					//sumU = L(j, k)*U(k, i);
					sumU += al[j][ku] * au[i][k];
				}

				//U(j, i) = A(j, i) - sumU;
				
					
					//L(i, j) = (A(i, j) - sumL) / U(j, j);
					if (i != j)
					{
						au[i][jl] = au[i][jl] - sumU;
						al[i][jl] = (al[i][jl] - sumL) / di[j];
					}
			}
			di[i] = di[i] - sumU;

		}
	}

	void FindY(form *y)
	{
		form *b = y;
		for (int i = 0; i < n; i++)
		{
			form2 sum = 0;
			for (int j = i-p, jl=0; j < i; j++, jl++)
			{
				if (j < 0) continue;
				//sum += L(i, j)*y[j];
				sum += al[i][jl] * y[j];
			}
			y[i] = b[i] - sum;
		}
	}

	void FindX(form *x)
	{
		form *y = x;
		for (int i = n - 1; i >= 0; i--)
		{
			
			x[i] = y[i] / di[i];
			/*if (i==n-1)
			{
				x[i] = 10;
			}*/
			for (int j = i - p, jl = 0; j < i; j++, jl++)
			{
				if (j < 0) continue;
				y[j] -= x[i] * au[i][jl];
			}
		}
	}
	/*form* Mult(const form* x, int n)
	{
		form* b = new form[n];
		for (size_t i = 0; i < n; i++)
		{
			b[i] = 0;
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = i-p, jl =0; j <=i; j++,jl++)
			{
				if (j < 0) continue;
				b[i] += al[i][jl] * x[j];
				b[j] += au[i][jl] * x[i];
			}
		}

		return b;
	}*/

	Matrix(int n)
	{
		this->n = n;
		p = n-1;
		m = 2 * p + 1;
		al = new form * [n];
		au = new form * [n];
		di = new form[n];

		for (int i = 0; i < n; i++)
		{
			al[i] = new form[p];
			au[i] = new form[p];

			for (int j = 0; j < p; j++)
			{
				al[i][j]=0;
				au[i][j]=0;
			}
			di[i]=0;
		}

	}
	Matrix(Matrix *M)
	{
		this->p = M->p;
		this->n = M->n;
		al = M->al;
		m = M->m;
		au = M->au;
		di = M->di;
	}
	Matrix(Matrix* M,int n, int p)
	{
		this->p = p;
		this->n = n;
		m = 2 * p + 1;
		al = new form * [n];
		au = new form * [n];
		di = new form[n];

		for (int i = 0; i < n; i++)
		{
			al[i] = new form[p];
			au[i] = new form[p];

			for (int j = 0; j < p; j++)
			{
				al[i][j]=M->al[i][j];
				au[i][j]=M->au[i][j];
			}
			di[i] = M->di[i];
		}
	}

	Matrix(std::ifstream &fl, std::ifstream &fu, std::ifstream &fd, int n, int k)
	{
		this->p = k;
		this->n = n;
		m = 2 * k + 1;
		al = new form*[n];
		au = new form*[n];
		di = new form[n];

		for (int i = 0; i < n; i++)
		{
			al[i] = new form[k];
			au[i] = new form[k];

			for (int j = 0; j < k; j++)
			{
				fl >> al[i][j];
				fu >> au[i][j];
			}
			fd >> di[i];
		}
	}


	~Matrix()
	{
		for (int i = 0; i < n; i++) {
			delete[]al[i];
			delete[]au[i];
		}

		delete[]al;
		delete[]au;
		delete[]di;
	}



	form & operator ()(const int &i, const int &j)
	{
		bool contains = false; //проверить
		if (abs(i - j) > p)
		{
			form tmp = 0;
			return tmp;
		}

		if (i == j) return di[i];
		if (j < i) return al[i][p - (i - j)];
		return au[j][p - (j - i)];
	}

	/*form operator ()(const int i, const int j) const
	{
		return 0;
	}*/

};

void Hilbert(Matrix &A, form* b)
{
	const double c1 = 1;
	for (int i = 0; i < A.n; i++)
	{
		for (int j = 0; j <=i; j++)
		{
			A(i, j) = c1 / (i + j + c1);
			A(j, i) = c1 / (i + j + c1);
		}
	}

	for (int i = 0; i < A.n; i++)
	{
		double sum = 0;
		for (size_t j = 0; j < A.n; j++)
		{
			sum += A(i, j) * (j + c1);
		}
		b[i] = sum;
	}
}

form** ForGauss( form* b, int n, int m)
{
	form **A = new form*[n];

	for (int i = 0; i < n; i++)
	{
		A[i] = new form[m];
		for (int j = 0; j < m; j++)
		{
			A[i][j] = 1.0 / (i + j + 1);
		}
	}

	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (size_t j = 0; j < n; j++)
		{
			sum += A[i][j] * (j + 1);
		}
		b[i] = sum;
	}
	return A;
}

void Gauss(form **A, form*b, int n, int m)
{
	for (size_t i = 0; i < n; i++)
	{
		form max = A[i][i];
		int l = i;
		for (int j = i; j < n; j++)
		{
			if (A[j][i] > max)
			{
				l = j;
				max = A[j][i];
			}
		}
		form *tmp = A[i];
		std::swap(A[i], A[l]);
		std::swap(b[i], b[l]);

		for (size_t j = i+1; j < n; j++)
		{
			form mult =A[j][i] / A[i][i];
			for (size_t k = i; k < n; k++)
			{
				if (k == i)
				{
					A[j][k] = 0;
				}
				else A[j][k] = mult*A[i][k];
			}
			b[j] -= b[i]*mult;
		}
	}

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			std::cout << A[i][j] << " ";
		}
		std::cout << "| " << b[i] <<'\n';
	}

	form *x = b;
	for (int i = n - 1; i >= 0; i--)
	{

		x[i] = b[i] / A[i][i];
		/*if (i==n-1)
		{
			x[i] = 10;
		}*/
		for (int j = 0; j < i; j++)
		{
			if (j < 0) continue;
			b[j] -= x[i] * A[i][j];
		}
	}
}
//class SLAE
//{
//public:
//	Matrix A;
//	Matrix L;
//	Matrix U;
//
//	form *b;
//	form *x;
//	form *y;
//
//	SLAE(Matrix& M, form *b):A(M), L(A), U(A)
//	{
//		this->b = b;
//		x = b;
//		y = b;
//	}
//
//	void Decompose()
//	{
//		for (int i = 0; i < A.n; i++)
//		{
//			for (int j = i - A.p; j <= i; j++)
//			{
//				if (j<0)
//				{
//					j = 0;
//				}
//				double sumL = 0;
//				double sumU = 0;
//
//				for (int k = i-A.p; k < j; k++)
//				{
//					if (k<0)
//					{
//						k = -1;
//						continue;
//					}
//					sumL = L(i, k)*U(k, j);
//					sumU = L(j, k)*U(k, i);
//				}
//
//				U(j, i) = A(j, i) - sumU;
//				
//				if (i == j) break;
//				L(i, j) = (A(i, j) - sumL) / U(j, j);
//			}
//		}
//	}
//
//	void findY()
//	{
//		for (int i = 0; i < A.n; i++)
//		{
//			double sum = 0;
//			for (int j = i-A.p; j<i; j++)
//			{
//				if (j < 0)
//				{
//					j = -1;
//					continue;
//				}
//				sum += L(i, j)*y[j];
//			}
//			y[i] = b[i] - sum;
//		}
//	}
//
//	void findX()
//	{
//		for (int i = A.n-1; i >=0 ; i--)
//		{
//			x[i] = y[i]/U(i,i);
//			for (int j =i-A.p; j<i; j++)
//			{
//				if (j<0)
//				{
//					j = -1;
//					continue;
//				}
//				y[j]-=x[i]*U(j,i);
//			}
//		}
//	}
//
//};




int main()
{
	int n, p;
	std::ifstream al;
	std::ifstream au;
	std::ifstream di;
	std::ifstream fb;
	std::ofstream tmp;
	tmp.open("out.txt");
	au.open("au.txt");
	al.open("al.txt");
	di.open("di.txt");
	fb.open("fb.txt");


	/*Matrix A = new Matrix(3);
	form *b = new form[3];
	Hilbert(A, b);*/


	////Hilbert Part
	form *b = new form[20];
	for (int k = 1; k < 20; k++)
	{
		
		Matrix A(k);
		Hilbert(A, b);
		A.Decompose();
		A.FindY(b);
		A.FindX(b);
		for (size_t i = 0; i < k; i++)
		{
			tmp << b[i]<<" "<< i+1 - b[i] << "\n";
		}

		tmp << "\n\n";
	}
	////----


	//fb >> n >> p;
	//

	/*Matrix A(al, au, di, n, p);
	form *x = new form[n];
	form* f = new form[n];
	for (size_t i = 0; i < n; i++)
	{
		fb >> f[i];
		x[i] = i + 1;
	}
	form* xk = new form[n];

	double multiplex = 1;
	for (int k = 0; k <=15; k++)
	{
		Matrix Ak(&A, n, p);
		for (size_t i = 0; i < n; i++)
		{
			xk[i] = f[i];
		}

		Ak(0, 0) += multiplex;
		for (size_t i = 0; i < n; i++)
		{
			xk[i] = f[i];
		}
		xk[0] += multiplex;
		Ak.Decompose();
		Ak.FindY(xk);
		Ak.FindX(xk);
		for (size_t i = 0; i < n; i++)
		{
			tmp << xk[i] << " " << x[i] - xk[i] << "\n";
		}
		multiplex /= 10;
	}*/

	//Gauss Part
	//form *b = new form[2];
	//form **A = ForGauss(b,2,2);
	//Gauss(A, b, 2, 2);
	//-------------
}
