#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace arma;
#define pi 3.14159265359

void compareArmaToAnalytical(int);
double offdiag(mat A, int* p, int* q, int n);
void bucklingBeam(int n);
void Jacobi_rotate(mat& A, mat& R, int k, int l, int n);
void harmonicOscillator(int n);
void twoElectrons(int n);

int main() {
    // Compile as 'c++ project.cpp -larmadillo'

    cout << "6 lowest eigenvalues of buckling beam problem:\n";
    int N = 200;   // Jacobi takes ages for N > 250
    compareArmaToAnalytical(N);
    bucklingBeam(N);

    cout << "\nHarmonic oscillator: \n";
    harmonicOscillator(N);

    cout << "\nTwo electrons, with interaction (unstable as a function of rho..): \n";
    twoElectrons(N);

    return 0;
}

void compareArmaToAnalytical(int N) {
    vec analyticalLambda(N);
    double h = 1. / N;
    double d = 2. / pow(h,2);
    double a = - pow(h,-2);

    for (int j = 1; j <= N; j++) {
        analyticalLambda(j-1) = d + 2*a*cos( j*pi / (N + 1) );
    }

    // Creating A
    vec aVec = zeros<vec>(N-1); aVec.fill(a);
    vec dVec = zeros<vec>(N); dVec.fill(d);
    mat A = diagmat(aVec,-1) + diagmat(dVec) + diagmat(aVec,1);
    vec lambdaArma; mat eigVecs;

    eig_sym(lambdaArma, eigVecs, A);    // Eigenvectors as columns in 'eigVecs'
    cout << "\nArmadillo solver - Analytical:\n";
    for (int i = 0; i < 6; i++) {
        cout << lambdaArma(i) << "           " << analyticalLambda(i) << endl;
    }
    // cout << eigVecs << endl;
}

void bucklingBeam(int n) {
    double h = 1. / n;
    double d = 2. / pow(h,2);
    double a = - pow(h,-2);
    vec aVec = zeros<vec>(n-1); aVec.fill(a);
    vec dVec = zeros<vec>(n); dVec.fill(d);

    mat A = diagmat(aVec,-1) + diagmat(dVec) + diagmat(aVec,1);
    mat R = eye<mat>(n,n);

    double tol = 1.e-10;
    int iteration = 0;
    int maxiter = pow(n,3);
    int p, q;
    double maxnondiag = offdiag(A, &p, &q, n);

    cout << "\nJacobi (N = " << n << "):\n";

    // Heavily influenced by the lecture notes
    while (maxnondiag > tol && iteration <= maxiter) {
        maxnondiag = offdiag(A, &p, &q, n);
        Jacobi_rotate(A, R, p, q, n);
        iteration++;
    }

    if (iteration > maxiter) {
        cout << "Nondiagonals still too large after max iteratoins\n";
    }

    // Eigenvalue stored on diagonal of A, and corresponding eigenvectors are stored in the columns(!) of R.
    // cout << A << endl;
    // cout << R << endl;

    vec lam = sort(diagvec(A));
    for (int i = 0; i < 6; i++) {
        cout << lam(i) << endl;
    }
}

void harmonicOscillator(int n) {
    double rho_max = 8.;
    rho_max = 4.;
    double h = rho_max / n;
    double d = 2. / pow(h,2);
    double a = - pow(h,-2);
    vec aVec = zeros<vec>(n-1); aVec.fill(a);
    vec dVec = zeros<vec>(n);
    for (int i = 0; i < n; i++) {
        dVec(i) = d + pow(i*h,2);
    }

    mat A = diagmat(aVec,-1) + diagmat(dVec) + diagmat(aVec,1);
    mat R = eye<mat>(n,n);

    // // Arma solver - allow for much higher n values
    // vec lambdaArma; mat eigVecs;
    // eig_sym(lambdaArma, eigVecs, A);
    // cout << "Eigenvalues:\n";
    // for (int i = 0; i < 5; i++) {
    //     cout << lambdaArma(i) << endl;
    // }
    // cout << ".\n" << ".\n" << ".\n" << endl;

    // Jacobi solver
    double tol = 1.e-10;
    int iteration = 0;
    int maxiter = pow(n,3);
    int p, q;
    double maxnondiag = offdiag(A, &p, &q, n);

    while (maxnondiag > tol && iteration <= maxiter) {
        maxnondiag = offdiag(A, &p, &q, n);
        Jacobi_rotate(A, R, p, q, n);
        iteration++;
    }

    vec lam = sort(diagvec(A));
    for (int i = 0; i < 6; i++) {
        cout << lam(i) << endl;
    }
}

void twoElectrons(int n) {
    double rho_max = 10.; // I'm not able to find a working value. The result vary regardless..
    double h = rho_max / n;
    double d = 2. / pow(h,2);
    double a = - pow(h,-2);
    vec aVec = zeros<vec>(n-1); aVec.fill(a);
    vec dVec = zeros<vec>(n);
    double omega_r[] = {.01, .5, 1., 5.};
    for (int i = 0; i < n; i++) {
        dVec(i) = d + pow(omega_r[0]*i*h,2) + pow(h*i,-1);
    }
    dVec(0) = 1e+15; // 1 / r, as r goes to zero

    mat A = diagmat(aVec,-1) + diagmat(dVec) + diagmat(aVec,1);
    mat R = eye<mat>(n,n);

    double tol = 1.e-10;
    int iteration = 0;
    int maxiter = pow(n,3);
    int p, q;
    double maxnondiag = offdiag(A, &p, &q, n);

    while (maxnondiag > tol && iteration <= maxiter) {
        maxnondiag = offdiag(A, &p, &q, n);
        Jacobi_rotate(A, R, p, q, n);
        iteration++;
    }

    vec lam = sort(diagvec(A));
    for (int i = 0; i < 6; i++) {
        cout << lam(i) << endl;
    }

    // ofstream outfile;
    // outfile.open("twoElectrons.txt");
    // for (int i = 0; i < n; i++) {
    //     // "rho, u(rho)" for ground state
    //     outfile << setprecision(8) << i*h << "," << R(i,0) << endl;
    // }
    // outfile.close();
}


// Assumes A is symmetric
void Jacobi_rotate(mat& A, mat& R, int k, int l, int n) {
    // COPY PASTE FROM FYS3150 GITHUB - did not have time to implement this myself. Don't give me any points here hehe
    // http://compphysics.github.io/ComputationalPhysics/doc/pub/eigvalues/html/eigvalues.html
    double s, c;
    if ( A(k,l) != 0.0 ) {
        double t, tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));

        if ( tau >= 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }

        c = 1./sqrt(1+t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
    A(l,k) = 0.0;  // same here
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
    //  And finally the new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);

        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}

double offdiag(mat A, int* p, int* q, int n) {
    // Copy pasta
   double max;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           {
              max = aij;  *p = i; *q = j;
           }
       }
   }
   return max;
}
