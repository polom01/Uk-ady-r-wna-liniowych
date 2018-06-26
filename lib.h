#pragma once
#include<iostream>
using namespace std;


class Matrix {
public:
	double **M;
	int N;
	Matrix(int N) {
		this->N = N;
		//create table
		M = new double *[N];
		for (int i = 0; i<N; i++)
			M[i] = new double[N];
		//insert value



	}
	void insert(int N, int a1, int a2, int a3) {

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (i == j) {
					M[i][j] = a1;
				}
				else if (i == j + 1 || i == j - 1) {
					M[i][j] = a2;
				}
				else if (i == j + 2 || i == j - 2) {
					M[i][j] = a3;
				}
				else {
					M[i][j] = 0;
				}

			}

		}
	}

	void print(int N) {
		for (int a = 0; a < N; a++) {
			for (int b = 0; b < N; b++) {
				cout << M[a][b] << " ";
			}
			cout << endl;
		}
	}

	int revert(int N) {

		for (int i = 0; i < N; i++) {
			if (M[i][i] == 0) {
				string wyjatek = "dzielenie przez zero!";
				throw wyjatek;
			}
			else {
				M[i][i] = 1 / M[i][i];
			}
		}


		return 1;
	}


	Matrix * operator + (Matrix const &oldMatrix) const {
		Matrix *retMatrix=new Matrix(this->N);
		
		for(int i=0;i<this->N;i++)
			for (int j = 0; j < this->N; j++) {
				retMatrix->M[i][j] = oldMatrix.M[i][j] + this->M[i][j];
		}
		return retMatrix;
	}

	Matrix * operator * (Matrix const &oldMatrix) const {
		Matrix *retMatrix = new Matrix(this->N);

		for (int i = 0; i < this->N; i++) {
			for (int j = 0; j < this->N; j++) {
				retMatrix->M[i][j] = 0;
				for (int k = 0; k < this->N; k++) {

					retMatrix->M[i][j] += this->M[i][k] * oldMatrix.M[k][j];
				}
			}
		}
		return retMatrix;
	}



	~Matrix()
	{
		//for (int i = 0; i < this->N; i++)
		//	delete M[i];
	//	delete M;
		
	}

};
