#include<iostream>
#include <fstream>
#include<math.h>
#include <time.h>
#include"lib.h"
#define LASTINDEXNUMBER 4
#define SECONDLASTINDEXNUMBER 2
#pragma warning(disable:4996)
using namespace std;

bool solve(int n, double ** L, double ** U, double * B, double * X);
double * createB(int N, int f) {
	double *B = new double[N];
	for (int i = 0; i < N; i++) {
		B[i] = sin(i*(f + 1));
	}
	return B;
}
void printB(int N, double *B) {
	for (int i = 0; i < N; i++) {
		cout << B[i];
	}
}
void jacobi(int N, double *B, Matrix matrix) {
	cout << "rozpoczyna metode jacobi" << endl;
	//M=L + D + U
	//creat table
	Matrix L(N);
	Matrix D(N);
	Matrix U(N);

	//insert value


	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				D.M[i][j] = matrix.M[i][j];
				L.M[i][j] = 0;
				U.M[i][j] = 0;
			}
			else if (i > j) {
				D.M[i][j] = 0;
				L.M[i][j] = 0;
				U.M[i][j] = matrix.M[i][j];
			}
			else if (i<j) {
				D.M[i][j] = 0;
				L.M[i][j] = matrix.M[i][j];
				U.M[i][j] = 0;
			}

		}
	}
	 D.revert(N);//D^-1
	//-D-1(L+U)
	 //xn+1=Mx+Nb
	// cout << "przed mnożeniem i dodawaniem" << endl;
	 Matrix *LU =L+U;
	// cout << "po mnożenia" << endl;
	// Matrix *multi = L * U;
	
	// M = -D ^-1(L + U)
	
	 Matrix M(N);
	 for (int i = 0; i<N; i++)
		 for (int j = 0; j<N; j++)
			 if (i == j)
				 M.M[i][j] = 0;
			 else
				 M.M[i][j] = -(LU->M[i][j] * D.M[i][i]);

	 //obliczenia 
	 double *X=new double[N];
	 double *X2=new double[N];
	 for (int i = 0; i < N; i++) {
		 X[i] = 0;//init
	 }
	 //xn+1=Mx+Nb

	 //petla
	 double r;
	 int count = 0;
	 double *residuum=new double[N];
	 do {
		 for (int i = 0; i < N; i++) {
			 X2[i] = D.M[i][i] * B[i];
			 for (int j = 0; j < N; j++) {

				 X2[i] += M.M[i][j] * X[j];
			 }
		 }
		 for (int i = 0; i < N; i++) {
			 X[i] = X2[i];
		 }

		 //r = A*x − b
		 for (int i = 0; i < N; i++) {
			 residuum[i] = 0;
			 for (int j = 0; j < N; j++) {
				 residuum[i] += matrix.M[i][j] * X[j];
			 }
			 residuum[i] -= B[i];
			
		 }
		 //norm
		 double norm = 0;
		 for (int i = 0; i < N; i++) {
			 norm += residuum[i] * residuum[i];
		 }
		r = sqrt(norm);
	
		count++;
		

	 } while (r>0.000000001);
	 for (int i = 0; i < N; i++) {
		// cout << "X(" << i << ") = " << X[i] << endl;
	 }
	
	 cout << "liczba powtórzeń : " << count << endl;;
	 delete LU;
	 //delete multi;
	 
}


void GaussaSeidla(int N, double *B, Matrix matrix){
	cout << "rozpoczony algorytm Gaussa–Seidla" << endl;
	//M=L + D + U
	//creat table
	Matrix L(N);
	Matrix D(N);
	Matrix U(N);

	//insert value


	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				D.M[i][j] = matrix.M[i][j];
				L.M[i][j] = 0;
				U.M[i][j] = 0;
			}
			else if (i > j) {
				D.M[i][j] = 0;
				L.M[i][j] = 0;
				U.M[i][j] = matrix.M[i][j];
			}
			else if (i<j) {
				D.M[i][j] = 0;
				L.M[i][j] = matrix.M[i][j];
				U.M[i][j] = 0;
			}

		}
	}
	D.revert(N);//D^-1
				//-D-1(L+U)
				//xn+1=Mx+Nb

	Matrix *LU = L + U;
//	Matrix *multi = L * U;

	// M = -D ^-1(L + U)

	Matrix M(N);
	for (int i = 0; i<N; i++)
		for (int j = 0; j<N; j++)
			if (i == j)
				M.M[i][j] = 0;
			else
				M.M[i][j] = -(LU->M[i][j] * D.M[i][i]);

	//obliczenia 
	double *X = new double[N];
	double *X2 = new double[N];
	for (int i = 0; i < N; i++) {
		X[i] = 0;//init
	}
	//xn+1=Mx+Nb

	//petla
	double r;
	int count = 0;
	double *residuum = new double[N];
	do {
		for (int i = 0; i < N; i++) {
			X2[i] = D.M[i][i] * B[i];
			for (int j = 0; j < N; j++) {
				
				X2[i] += M.M[i][j] * X[j];
			}
			X[i] = X2[i];
		}

		//r = A*x − b
		for (int i = 0; i < N; i++) {
			residuum[i] = 0;
			for (int j = 0; j < N; j++) {
				residuum[i] += matrix.M[i][j] * X[j];
			}
			residuum[i] -= B[i];

		}
		//norm
		double norm = 0;
		for (int i = 0; i < N; i++) {
			norm += residuum[i] * residuum[i];
		}
		r = sqrt(norm);

		count++;


	} while (r>0.000000001);
	for (int i = 0; i < N; i++) {
		//cout << "X(" << i << ") = " << X[i] << endl;
	}

	cout << "liczba powtórzeń : " << count << endl;;

	delete LU;



}
Matrix* addTwoMatrix(Matrix *left,Matrix *right) {

	Matrix *retMatrix= new Matrix(right->N);
	for (int i = 0; i<right->N; i++)
		for (int j = 0; j < right->N; j++) {
			retMatrix->M[i][j] = left->M[i][j] + right->M[i][j];
		}
	return retMatrix;
}


void zadB() {

	cout << "zadanie B" << endl;

	int e = 3;
	int f = 5;
	int N = 9 * LASTINDEXNUMBER*SECONDLASTINDEXNUMBER;
	int a1 = 5 + e;
	int	a2 = -1;
	int	a3 = -1;
	Matrix jacobiMatrix(N);
	jacobiMatrix.insert(N, a1, a2, a3);
	double *B = createB(N, f);
	jacobi(N, B, jacobiMatrix);
	cout << endl << endl;

	cout << N << endl;


	Matrix GaussMatrix(N);
	GaussMatrix.insert(N, a1, a2, a3);
	GaussaSeidla(N, B, GaussMatrix);
	return;

}

void factorization(Matrix *matrix,double * X,double *b) {

	Matrix I(matrix->N);
	for (int a = 0; a < matrix->N; a++) {
		for (int b = 0; b < matrix->N; b++) {
			if (a == b)
				I.M[a][b] = 1;
			else
				I.M[a][b] = 0;
		}
	}
	for (int k = 1; k < matrix->N-1; k++) {
		for (int j = k + 1; j <matrix->N; j++) {
			I.M[j][k] = matrix->M[j][k] / matrix->M[k][k];
			for (int z = k; z < matrix->N; z++) {
				matrix->M[j][z] = matrix->M[j][z] - matrix->M[k][z]*I.M[j][k];
			}
		
		
		}
	}
	solve(matrix->N,I.M, matrix->M, b, X);
	for (int a = 0; a < matrix->N; a++) {
		//cout << "X(" << a << ") = " << X[a] << endl;
	}
	
	return;
}

bool solve(int n, double ** L,double ** U, double * B, double * X)
{
	double s;

	X[0] = B[0];

	for (int i = 1; i < n; i++)
	{
		s = 0;
		for (int j = 0; j < i; j++) s += L[i][j] * X[j];
		X[i] = B[i] - s;
	}


	X[n - 1] /= U[n - 1][n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		s = 0;
		for (int j = i + 1; j < n; j++) s += U[i][j] * X[j];
		X[i] = (X[i] - s) / U[i][i];
	}

	return true;
}


void gauss(Matrix *matrix,double *B) {
	cout << endl << "gauss" << endl;
	int N = matrix->N;
	//matrix->print(N);
	for (int x = 0; x < N - 1; x++) {//N-1 ostatniego nie mozna zrobic
		int row = x + 1;


		for (int z = x; z < N-1; z++) {
			double howMuch = 0;
			if (matrix->M[x][x] != 0)
				howMuch = matrix->M[z+1][x] / matrix->M[x][x];
			for (int y = x; y < N; y++) {
				matrix->M[z + 1][y] -= howMuch * matrix->M[x][y];//dodac na odejmowanie B
			}
			B[z + 1] -= howMuch*B[z];
		}
	}
	//matrix->print(N);


	//obliczenia podstawiienie wstecz
	double *X = new double[N];
	X[N-1] = B[N-1]/ matrix->M[N-1][N-1];
	//cout << endl << B[N - 1];
	//cout << endl << matrix->M[N - 1][N - 1];
	//cout << endl << endl << X[N - 1] << endl << endl;
	for (int i = N - 2; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < N; j++)sum += matrix->M[i][j] * X[j];
		X[i] = (B[i] - sum) / matrix->M[i][i];
	}

	for (int i = 0; i < N; i++) {
		//cout << "X(" << i << ") = " << X[i] << endl;
	}
}



void zadC() {

	cout << endl<<"zadanie C" << endl;
	int e = 3;
	int f = 5;
	int N = 9 * LASTINDEXNUMBER*SECONDLASTINDEXNUMBER;
	int a1 = 3;
	int	a2 = -1;
	int	a3 = -1;
	Matrix jacobiMatrix(N);
	jacobiMatrix.insert(N, a1, a2, a3);
	double *B = createB(N, f);
	jacobi(N, B, jacobiMatrix);
	cout << endl << endl;

	//cout << N << endl;
	//printB(N, B);

	Matrix GaussMatrix(N);
	GaussMatrix.insert(N, a1, a2, a3);
	GaussaSeidla(N, B, GaussMatrix);


	return;
}

void zadD() {

	cout << endl << "zadanie D" << endl;
	int e = 3;
	int f = 5;
	int N = 9 * LASTINDEXNUMBER*SECONDLASTINDEXNUMBER;
	int a1 = 5+e;
	int	a2 = -1;
	int	a3 = -1;
	Matrix LuMatrix(N);
	LuMatrix.insert(N, a1, a2, a3);
	double *B = createB(N, f);
	//Ax=b
	double *X = new double[LuMatrix.N];
	factorization(&LuMatrix,X,B);


	Matrix GaussMatrix(N);
	GaussMatrix.insert(N, a1, a2, a3);
	double *B2 = createB(N, f);
	gauss(&GaussMatrix,B2);

}

void zadE() {
	int e = 3;
	int f = 5;
	int N[5] = { 100, 500, 1000, 2000, 3000 };
	int a1 = 5 + e;
	int	a2 = -1;
	int	a3 = -1;
	ofstream myFile;
	myFile.open("data.txt");
	
	

	//jacobi
	cout << endl << "jacobi" << endl;
	for (int i = 0; i < 5; i++) {
		Matrix matrix(N[i]);
		matrix.insert(N[i], a1, a2, a3);
		double *B = createB(N[i], f);

		int start = clock();
		jacobi(matrix.N, B,matrix);
		double time = clock() - start;
		cout << "dla N = " << N[i] << " czas : " << time << endl;
		myFile << "jacobi" << endl;
		myFile << time<<endl;

	}
	
	//gaussaseidla
	cout << endl << "GaussaSeidla" << endl;
	for (int i = 0; i < 5; i++) {
		Matrix matrix(N[i]);
		matrix.insert(N[i], a1, a2, a3);
		double *B = createB(N[i], f);

		int start = clock();
		GaussaSeidla(matrix.N, B, matrix);
		double time = clock() - start;
		cout << "dla N = " << N[i] << " czas : " << time << endl;
		myFile << "GaussaSeidla" << endl;
		myFile << time << endl;


	}

	//faktoryzacja LU

	cout << endl << "faktoryzacja LU" << endl;
	for (int i = 0; i < 5; i++) {
		Matrix matrix(N[i]);
		matrix.insert(N[i], a1, a2, a3);
		double *B = createB(N[i], f);
		double *X = new double[matrix.N];
		int start = clock();
		factorization(&matrix, X, B);
		double time = clock() - start;
		cout << "dla N = " << N[i] << " czas : " << time << endl;
		myFile << "faktoryzacja LU" << endl;
		myFile << time << endl;

	}


	//gauss
	cout << endl << "gauss" << endl;
	for (int i = 0; i < 5; i++) {
		Matrix matrix(N[i]);
		matrix.insert(N[i], a1, a2, a3);
		double *B = createB(N[i], f);

		int start = clock();
		gauss(&matrix, B);
		double time = clock() - start;
		cout << "dla N = " << N[i] << " czas : " << time << endl;
		myFile << "gauss" << endl;
		myFile << time << endl;


	}
	
	myFile.close();
}

int main() {
	

	
	
	zadB();
	zadC();
	zadD();
	zadE();
	
	int e;
	cin >> e;
	return 0;
}