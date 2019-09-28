﻿// NM1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <math.h>
#include <cstdio>
#include <fstream>
#include <iostream>

typedef float form;

void Test(std::ifstream &fi)
{

}

class Matrix
{
public:
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
			for (int j = i-p, jl=0; jl <= p&&j<=i; j++, jl++)
			{
				if (j < 0) continue;
				double sumL = 0;
				double sumU = 0;

				for (int k = 0,ku=i-j; k < jl; k++, ku++)
				{
					//sumL = L(i, k)*U(k, j);
					sumL += al[i][k] * au[j][ku];
					//sumU = L(j, k)*U(k, i);
					sumU += al[j][ku] * au[i][k];
				}

				//U(j, i) = A(j, i) - sumU;
				
				if (i == j)
				{
					di[j] = di[j] - sumU;
					break;
				}
				else
				{
					au[i][jl] = au[i][jl] - sumU;
				}
				//L(i, j) = (A(i, j) - sumL) / U(j, j);
				al[i][jl] = (al[i][jl] - sumL) / di[j];
			}
		}
	}

	void FindY(form *y)
	{
		form *b = y;
		for (int i = 0; i < n; i++)
		{
			double sum = 0;
			for (int j = i-p, jl=0; j < i; j++)
			{
				if (j < 0) continue;
				//sum += L(i, j)*y[j];
				sum += al[i][jl] * b[j];
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
			for (int j = i - p, jl = 0; j < i; j++, jl++)
			{
				if (j < 0) continue;
				y[j] -= x[i] * au[i][jl];
			}
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

	/*~Matrix()
	{
		for (int i = 0; i < n; i++)
		{
			delete al[i];
			delete au[i];
		}
		delete di;
		delete al;
		delete au;
	}*/


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

class SLAE
{
public:
	Matrix A;
	Matrix L;
	Matrix U;

	form *b;
	form *x;
	form *y;

	SLAE(Matrix& M, form *b):A(M), L(A), U(A)
	{
		this->b = b;
		x = b;
		y = b;
	}

	void Decompose()
	{
		for (int i = 0; i < A.n; i++)
		{
			for (int j = i - A.p; j <= i; j++)
			{
				if (j<0)
				{
					j = 0;
				}
				double sumL = 0;
				double sumU = 0;

				for (int k = i-A.p; k < j; k++)
				{
					if (k<0)
					{
						k = -1;
						continue;
					}
					sumL = L(i, k)*U(k, j);
					sumU = L(j, k)*U(k, i);
				}

				U(j, i) = A(j, i) - sumU;
				
				if (i == j) break;
				L(i, j) = (A(i, j) - sumL) / U(j, j);
			}
		}
	}

	void findY()
	{
		for (int i = 0; i < A.n; i++)
		{
			double sum = 0;
			for (int j = i-A.p; j<i; j++)
			{
				if (j < 0)
				{
					j = -1;
					continue;
				}
				sum += L(i, j)*y[j];
			}
			y[i] = b[i] - sum;
		}
	}

	void findX()
	{
		for (int i = A.n-1; i >=0 ; i--)
		{
			x[i] = y[i]/U(i,i);
			for (int j =i-A.p; j<i; j++)
			{
				if (j<0)
				{
					j = -1;
					continue;
				}
				y[j]-=x[i]*U(j,i);
			}
		}
	}

};




int main()
{
	int n, k;
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

	fb >> n >> k;
	form **test = new form*[n];
	Matrix A(al,au,di,n,k);
	Matrix B(A);
	form *b = new form[n];
	for (size_t i = 0; i < n; i++)
	{
		fb >> b[i];
	}

	A.Decompose();
	A.FindY(b);

	//SLAE S(A,b);
	//S.Decompose();
	//S.findY();
	//S.findX();
	//int t = 0;
	
}
