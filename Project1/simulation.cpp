#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <armadillo>

using namespace std;
using namespace arma;

double f(double);
double u_analytic(double);
void generalized_tridiagonal_solver(int n);
void write_to_file(string filenam, double* x, double* u, double* v, int n);
double identical_diagonals_solver(int n, bool write);
double find_max(double* lst, int n);
void general_matrix_solver(int n);


int main(int argc, char* argv[]) {

    if (argc > 1) {
        // User-defined n, and write results to file.
        int n;  // Number of interior points
        double tmp = atof(argv[1]); // This makes sure we can read the number in scientific notation.
        n = (int) tmp; // Needs to be int though.
        generalized_tridiagonal_solver(n);
        bool write = true; // writes solution to file
        identical_diagonals_solver(n, write); // Set 'true' for writing results to file.
        general_matrix_solver(n); // out of memory at n = 1e5..
    } else {
        // Solves for different values of n and write errors to file.
        // I only use the fastest solver for this
        ofstream outfile;
        outfile.open("errors.txt");

        int n_lst[] = {10,100,1000,10000,100000,1000000,10000000};
        double log_h, log_eps;
        for (int n : n_lst) {
            log_h = log10(1. / (n + 1.));

            bool write = false; // don't write solution to file
            log_eps = identical_diagonals_solver(n, write);

            // Write errors to file. Format: "log10(h),log10(eps)\n"
            outfile << setprecision(8) << log_h << "," << log_eps << endl;
        }
        outfile.close();
    }

    // This got a bit messy hehe, classes next time :))
    return 0;
}

double f(double x) {
    return 100. * exp(-10.*x);
}

double u_analytic(double x) {
    return 1 - (1 - exp(-10.))*x - exp(-10.*x);
}

void generalized_tridiagonal_solver(int n) {
    double h = 1. / (n + 1);    // Step length

    // Allocate memory at run time.
    double* a = new double[n+2];
    double* b = new double[n+2];
    double* c = new double[n+2];
    double* v = new double[n+2];    // Numeric solution
    double* u = new double[n+2];    // Analytical solution
    double* x = new double[n+2];
    double* b_tilde = new double[n+2];

    // Initialize vectors. Using same indexing for all
    v[0] = 0; v[n+1] = 0;   // Boundary condition
    for (int i = 0; i <= n+1; i++) {
        x[i] = h*i;
        u[i] = u_analytic(x[i]);    // Analytical
        if (i <= n && i >= 1) {
            b[i] = 2;
            if (i < n) {
                a[i] = -1;
                c[i] = -1;
            }
            b_tilde[i] = pow(h,2.)*f(x[i]);
        }
    }

    // Gaussian elimination
    clock_t start, finish;
    start = clock();
    for (int i = 2; i <= n; i++) {
        b[i] -= c[i-1] * (a[i-1]/b[i-1]);               // 3*(n-1) flops
        b_tilde[i] -= b_tilde[i-1] * (a[i-1]/b[i-1]);     // 3*(n-1) flops
    }
    v[n] = b_tilde[n] / b[n];                           // 1 flop
    for (int i = (n-1); i >= 1; i--) {
        v[i] = (b_tilde[i] - c[i]*v[i+1]) / b[i];       // 3*(n-1) flops
    }
    finish = clock();
    double time_used = (double)(finish - start)/CLOCKS_PER_SEC * 1000.; // [ms]
    cout << "General tri. dig.: " << time_used << " ms" << endl;

    write_to_file("simulation1.txt", x, u, v, n);

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] u;
    delete[] v;
    delete[] x;
    delete[] b_tilde;
}

double identical_diagonals_solver(int n, bool write) {
    // This method also calculates the error

    double h = 1. / (n + 1);    // Step length
    double a = -1;
    // double b = 2;
    double c = -1;

    // Allocate memory at run time.
    double* b = new double[n+2];
    double* x = new double[n+2];
    double* u = new double[n+2];    // Analytical solution
    double* v = new double[n+2];    // Numeric solution
    double* b_tilde = new double[n+2];

    // Initialize vectors. Using same indexing for all
    v[0] = 0; v[n+1] = 0;   // Boundary condition
    for (int i = 0; i <= n+1; i++) {
        x[i] = h*i;
        u[i] = u_analytic(x[i]);    // Analytical
        if (i <= n && i >= 1) {
            b[i] = 2;
            b_tilde[i] = pow(h,2.)*f(x[i]);
        }
    }

    // Gaussian elimination
    double ac = a * c;
    clock_t start, finish;
    start = clock();
    for (int i = 2; i <= n; i++) {
        b[i] -= ac / b[i-1];
        b_tilde[i] -= b_tilde[i-1] * (a/b[i-1]);
    }
    v[n] = b_tilde[n] / b[n];
    for (int i = (n-1); i >= 1; i--) {
        v[i] = (b_tilde[i] - c*v[i+1]) / b[i];
    }


    finish = clock();
    double time_used = (double)(finish - start)/CLOCKS_PER_SEC * 1000.; // [ms]
    cout << "Identical diagonals: " << time_used << " ms" << endl;

    if (write) {
        write_to_file("simulation2.txt", x, u, v, n);
    }

    // Rel error for interior points
    double* rel_error = new double[n];
    for (int i = 1; i <= n; i++) {
        rel_error[i-1] = log10(fabs((v[i] - u[i]) / u[i]));
    }
    double max = find_max(rel_error, n);

    delete[] x;
    delete[] u;
    delete[] v;
    delete[] b_tilde;

    return max;
}

void general_matrix_solver(int n) {
    vec a = zeros<vec>(n-1); a.fill(-1);
    vec b = zeros<vec>(n); b.fill(2);
    vec c = zeros<vec>(n-1); c.fill(-1);
    mat A = diagmat(a,-1) + diagmat(b) + diagmat(c,1);

    vec b_tilde = zeros<vec>(n);
    vec x = zeros<vec>(n);
    double h = 1. / (n + 1);
    for (int i = 0; i < n; i++) {
        x(i) = (i + 1) * h;
        b_tilde(i) = pow(h,2.) * f(x(i));
    }
    // Solver
    clock_t start, finish;
    start = clock();
    vec v = solve(A, b_tilde);
    finish = clock();

    double time_used = (double)(finish - start) / CLOCKS_PER_SEC * 1000.; // ms
    cout << "LU: " << time_used << " ms.\n";

    v.save("v.txt", raw_ascii);
    x.save("x.txt", raw_ascii);
}

void write_to_file(string filename, double* x, double* u, double* v, int n) {
    ofstream outfile;
    outfile.open(filename);

    for (int i = 0; i <= n+1; i++) {
        // Format: "x,u,v\n"
        outfile << setprecision(8) << x[i] << "," <<
        u[i] << "," << v[i] << endl;
    }
    outfile.close();
}

double find_max(double* lst, int n) {
    double current_max = lst[0];
    for (int i = 1; i < n; i++) {
        if (lst[i] > current_max) {
            current_max = lst[i];
        }
    }
    return current_max;
}
