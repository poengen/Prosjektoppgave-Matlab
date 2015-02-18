//C++ code for solving poisson 2D

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char** argv){
	
	// THE GRID
	int n;
	float h;
	n = 3;
	h = 1/float(n+1);

	vec x = zeros<vec>(n+2);
	for (int i = 0; i<n+2; i++){
		x(i) = i*h;
	}
	
	//cout << h << endl;
	
	// MATRICES
	mat A = zeros<mat>(n*n,n*n);

	for (int i = 0; i < n*n; i++){
		A(i,i) = 4;
	}
	
	for (int i = 0; i < n*n-n; i++) {
		A(i,n+i) = -1;
		A(n+i,i) = -1;
	}

	for (int i = 0; i < n*n-1; i++) {
		A(i,i+1) = -1;
		A(i+1,i) = -1;
	}

	for (int i = 0; i < n-1; i++) {
		A((i+1)*n-1,(i+1)*n) = 0;
		A((i+1)*n,(i+1)*n-1) = 0;
	}

	cout << A << endl;
	
	//	solve linear system
	
	
	
	
	
	// compose solution matrix
	
	
	
	
	// save solution
	//Usol.save("A.dat", raw_ascii);
	
}
