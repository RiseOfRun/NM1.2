// NM1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
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
	int k;

	Matrix(Matrix *M)
	{
		this->k = M->k;
		this->n = M->n;
		al = M->al;
		m = M->m;
		au = M->au;
		di = M->di;
	}

	Matrix(std::ifstream &fl, std::ifstream &fu, std::ifstream &fd, int n, int k)
	{
		this->k = k;
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
		for (int i = 0; i < n; i++)
		{
			delete al[i];
			delete au[i];
		}
		delete di;
		delete al;
		delete au;
	}


	form & operator ()(const int &i, const int &j)
	{
		bool contains = false; //проверить
		if (abs(i - j) > k)
		{
			form tmp = 0;
			return tmp;
		}

		if (i == j) return di[i];
		if (j < i) return al[i][k - (i - j)];
		return au[j][k - (j - i)];
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
			for (int j = i - A.k; j <= i; j++)
			{
				if (j<0)
				{
					j = 0;
				}
				double sumL = 0;
				double sumU = 0;

				for (int k = i-A.k; k < j; k++)
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
			for (int j = i-A.k; j<i; j++)
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
			for (int j =i-A.k; j<i; j++)
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

	Matrix A(al,au,di,n,k);
	Matrix B(A);
	form *b = new form[n];
	for (size_t i = 0; i < n; i++)
	{
		fb >> b[i];
	}

	SLAE S(A,b);
	S.Decompose();
	S.findY();
	S.findX();
	int t = 0;
	
}
