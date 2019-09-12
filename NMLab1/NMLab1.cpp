// NM1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <math.h>
#include <cstdio>
#include <fstream>

typedef float form;


class System
{
	Matrix A;
	Matrix L;
	Matrix U;

	form *b;
	form *x;
	form *y;

	System(Matrix M, form *b):A(M), L(A), U(A)
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
				double sumL = 0;
				double sumU = 0;

				for (int k = i - A.k; k < i; k++)
				{
					sumU += L(j, k)*U(k, i);
				}
				for (int k = i - A.k; k < j; k++)
				{
					sumL += L(i, k)*U(k, j);
				}

				U(j, i) = A(j, i) - sumU;
				
				if (i == j) break;
				L(i, j) = (A(i, j) - sumL) / U(j, j);
			}
		}
	}

	void findY()
	{
		for (int i = 0; i < n; i++)
		{

		}
	}
};

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

	Matrix(std::ifstream fl, std::ifstream fu, std::ifstream fd, int n, int k)
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


	form & operator ()(const int &i, const int &j) const
	{
		bool contains = false; //проверить
		if (abs(i-j) > k)
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


int main()
{
	M(2, 2) = 4;
	Matrix L(M);
	Matrix U(M);
	M(2, 2) = 5;
	for (int i = 0; i < M.n; i++)
	{
		form sum = 0;
		for (int k = abs(i - (M.k + 1)); k < i; k++)
		{
			sum += L(i, k)*U(k, i);
		}
		U(i, i) = M(i, i) - sum;
		for (int j = i + 1; j < i + M.k&&j < M.n; j++)
		{
			sum = 0;
			form sum2 = 0;
			for (int k = i - M.k + 1; k < i; k++)
			{
				sum += L(i, k)*U(k, j);
			}
			U(i, j) = M(i, j) - sum;

			sum = 0;
			for (int k = i - M.k + 1; k < i; k++)
			{
				sum += L(j, k)*U(k, i);
			}
			L(j, i) = (M(j, i) - sum) / U(i, i);
		}
	}
}
