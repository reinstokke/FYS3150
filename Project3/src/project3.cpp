#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>
#include <time.h>
#include <fstream>
#include <iomanip>
// #include "ODESolver.h"

using namespace std;
using namespace arma;

class ODESolver {
    // I was not able to make the class in a seperate file. Resulted in "undefined reference [...]" error. Was not able to solve it.

    private:
        int N;                  // Number of steps
        double h;               // Step length. h = dt [yrs] in this case
        double beta;            // Parameter for the force: 1 / r^beta
        double G;               // Gravitational constant
        bool inertial_sun;      // True if the Sun is not accelerated. We allow it to move though. (v=v0=const)
        int dim;                // Spatial dimension of the problem
        bool mercury_GR;

    public:
        ODESolver(int dim=2, int N=1, double h=0, double beta=2., bool merc_GR=false) {
            this->N = N;
            this->h = h;
            G = 4. * pow(M_PI,2);   // [M] = solar masses
            inertial_sun = false; // This is false by default, but can be altered in a setter method.
            this->beta = beta;
            this->dim = dim;
            mercury_GR = merc_GR;
        }

    private:
        class Body {
            // Nested class for the bodies in the system. I don't feel like the user of the ODESolver should have to create objects for all the celestial bodies to use it. It was at least an attempt to hide the complexity from the user.

            public:
                double mass;
                double mass_solar; // Just the planet's mass in units of solar masses!
                string name;
                rowvec r0, v0;
                mat r, v;          // Stores the positions and velocities for all times.

                Body(double mass, string name, rowvec r0, rowvec v0, int N, int dim) {
                    this->name = name;
                    this->mass = mass; mass_solar = mass / (2.E+30);
                    if (name == "Sun") {
                        // Just to be exact
                        mass_solar = 1.;
                    }
                    this->r0 = r0; this->v0 = v0;
                    r = zeros<mat>(N,dim);
                    v = zeros<mat>(N,dim);
                    r.row(0) = r0;
                    v.row(0) = v0;
                }

                void setMass(double mass) {
                    this->mass = mass;
                    mass_solar = mass / 1.989E+30;
                }
        };

        vector<Body*> bodies;   // Stores all the body objects (pointers) in the solver

    public:
        void add_body(string name, double mass, rowvec r0, rowvec v0) {
            bodies.push_back(new Body(mass, name, r0, v0, N, dim));
        }

        void set_resolution(int N, double h) {
            // Unfortunate way to implement it. But if I want to create v and r matrices, I need to know the size at instantiation of the object. But I want to be able to change the resolution later too.

            this->N = N;
            this->h = h;
            for (Body* body: bodies) {
                body->r.reset();
                body->v.reset();
                body->r = zeros<mat>(N,dim);
                body->v = zeros<mat>(N,dim);
                body->r.row(0) = body->r0;
                body->v.row(0) = body->v0;
            }
        }

        void initialize_body(string name, rowvec r0, rowvec v0) {
            // If I want to change the initial values of a body after adding it to the solver.
            for (Body* body : bodies) {
                if (body->name == name) {
                    body->r0 = r0;
                    body->v0 = v0;
                    body->r.row(0) = body->r0;
                    body->v.row(0) = body->v0;
                }
            }
        }

        void set_body_mass(string name, double mass) {
            for (Body* body : bodies) {
                if (body->name == name) {
                    body->setMass(mass);
                }
            }
        }

        void euler() {
            // Note: Not Euler Cromer
            Body* current = NULL;
            Body* other = NULL;
            rowvec a = zeros<rowvec>(dim);  // Holds the acceleration for a given body and a given iteration

            for (int i = 0; i < (N - 1); i++) {
                for (int j = 0; j < bodies.size(); j++) {
                    current = bodies[j];

                    // Calculate acceleration
                    a.zeros(); // Resets betweem each iteration
                    if (!inertial_sun || current->name != "Sun") {
                    // Acceleration remains zero for inertial sun
                        for (int k = 0; k < bodies.size(); k++) {
                            other = bodies[k];
                            if (other->name != current->name) {
                                a += acceleration_gravity(current->r.row(i), other->r.row(i), other->mass_solar);
                            }
                        }
                    }

                    // I actually update r_sun even if the Sun is not accelerated, as it would move in a straight line if the initial velocity differ from zero. This is unnecessary if the initial velocity is zero, so in that case this solver performs many unnecessary operations. I do the same for the Verlet solver.
                    current->v.row(i+1) = current->v.row(i) + h * a;
                    current->r.row(i+1) = current->r.row(i) + h * current->v.row(i);
                }
            }
        }

        void velocity_verlet() {

            mat A1 = zeros<mat>(bodies.size(),dim); // Stores acceleration for all objects at time i
            mat A2 = zeros<mat>(bodies.size(),dim); // same for i+1

            Body* current;
            Body* other;

            // Calculate initial acceleration for all bodies
            for (int j = 0; j < bodies.size(); j++) {
                current = bodies[j];
                if (current->name == "Sun" && inertial_sun) {
                    // Acceleration remains zero
                    continue;
                }
                if (current->name == "Mercury" && mercury_GR) {
                    A1.row(j) = acceleration_mercury_GR(current->r.row(0), current->v.row(0));
                    continue;
                }
                for (int k = 0; k < bodies.size(); k++) {
                    other = bodies[k];
                    if (current->name != other->name) {
                        A1.row(j) += acceleration_gravity(current->r.row(0), other->r.row(0), other->mass_solar);
                    }
                }
            }

            // Moving some operations out of the loop
            double hh_half = pow(h,2)/2.;
            double h_half = h / 2.;

            // Starting 'velocity Verlet':
            for (int i = 0; i < (N - 1); i++) {

                // Calculates r(i+1) for all bodies
                for (int j = 0; j < bodies.size(); j++) {
                    current = bodies[j];
                    current->r.row(i+1) = current->r.row(i) + h*current->v.row(i) + hh_half * A1.row(j);
                }

                // Calculate A2
                A2.zeros(); // This made me debug for an hour...
                for (int j = 0; j < bodies.size(); j++) {
                    current = bodies[j];
                    if (current->name == "Sun" && inertial_sun) {
                        continue;
                    }
                    if (current->name == "Mercury" && mercury_GR) {
                        // Added this in the end. I hope everything else still works..
                        A2.row(j) = acceleration_mercury_GR(current->r.row(i+1), current->v.row(i+1));
                        continue;
                    }
                    for (int k = 0; k < bodies.size(); k++) {
                        other = bodies[k];
                        if (current->name != other->name) {
                            A2.row(j) += acceleration_gravity(current->r.row(i+1), other->r.row(i+1), other->mass_solar);
                        }
                    }
                }

                // Calculates v(i+1) for all bodies
                for (int j = 0; j < bodies.size(); j++) {
                    current = bodies[j];
                    current->v.row(i+1) = current->v.row(i) + h_half * (A2.row(j) + A1.row(j));
                }

                A1 = A2; // Faster to swap these in each iteration than to calculate the acceleration twice in every iteration.
            }
        }

        void write_to_file(string name="none", string filename="none") {
            if (name == "none") {
                for (int i = 0; i < bodies.size(); i++) {
                    string file = bodies[i]->name + ".dat";
                    bodies[i]->r.save(file, raw_ascii);
                }
                cout << "Data files written out as '<planet_name>.dat'.\n";
            }
            else {
                for (Body* body : bodies) {
                    if (body->name == name) {
                        body->r.save(filename, raw_ascii);
                    }
                }
                string msg = "Positions for " + name + " is written to '" + filename + "'.";
                cout << msg << endl;
            }
        }

        void setBeta(double b) {
            beta = b;
        }

        void set_inertial_sun(bool flag) {
            inertial_sun = flag;
        }

        double check_stability_circular(bool print=false) {
            // Test if r, v**2 and r x v are all constant.
            // Returning log10(max(abs(r - r0)))

            if (print) {
                double dt_seconds = h * 365.*24*60*60;
                double T = h*N;
                cout << "Stability check for circular orbit:\n";
                cout << "Simulation period = " << T << " yrs, " << to_string(N) << " steps.\n";
                string msgg = "dt = " + to_string(dt_seconds) + " seconds.\n";
                cout << msgg;
            }

            // Just finding the earth pointer
            Body* earth;
            for (Body* body : bodies) {
                if (body->name == "Earth") {
                    earth = body;
                }
            }

            double r_init = norm(earth->r0,2);
            double v2_init = pow(norm(earth->v0,2), 2);
            double ang_init = abs(earth->r0(0)*earth->v0(1) - earth->r0(1)*earth->v0(0));

            double r_error = 0;
            double v2_error = 0;
            double ang_error = 0;

            double tmp1, tmp2, tmp3;
            for (int i = 1; i < N; i++) {
                tmp1 = abs( norm(earth->r.row(i),2) - r_init );
                tmp2 = abs( pow(norm(earth->v.row(i),2), 2) - v2_init );
                tmp3 = abs( abs(earth->r.row(i)(0)*earth->v.row(i)(1) - earth->r.row(i)(1)*earth->v.row(i)(0)) - ang_init);

                // Finding biggest errors
                if (tmp1 > r_error) {
                    r_error = tmp1;
                }
                if (tmp2 > v2_error) {
                    v2_error = tmp2;
                }
                if (tmp3 > ang_error) {
                    ang_error = tmp3;
                }
            }

            // Expressed as relative errors
            r_error = log10( r_error / r_init );
            v2_error = log10( v2_error / v2_init );
            ang_error = log10( ang_error / ang_init );

            if (print) {
                cout << "log10 of relative error of position = " << r_error << endl;
                cout << "log10 of relative error of kinetic energy = " << v2_error << endl;
                cout << "log10 of relative error of angular momentum = " << ang_error << endl;
            }

            return r_error;
        }

        bool check_if_enclosed() {
            Body* earth;
            for (Body* body : bodies) {
                if (body->name == "Earth") {
                    earth = body;
                }
            }

            double eps = 5e-2; // Not very strict, because if we do reach escape velocity the planet will not even come close to the starting point again.
            bool bound = false;
            bool flag = false;
            for (int i = 1000; i < N; i++) {
                if (earth->r.row(i)(1) > eps && !flag) {
                    flag = true;
                }
                if (abs(earth->r0(0) - earth->r.row(i)(0)) < eps
                && abs(earth->r0(1) - earth->r.row(i)(1)) < eps
                && flag) {
                    // Back at same position - closed elliptic orbit
                    bound = true;
                    cout << i << endl;
                    break;
                }
            }
            return bound;
        }

        void set_mercury_GR(bool flag) {
            mercury_GR = flag;
        }

        void calculate_precession() {
            // I'll only calculate the theta for the last time it passes the perihelion, and then return theta_p_100 - theta_p_0

            Body* mercury;
            for (Body* body : bodies) {
                if (body->name == "Mercury") {
                    mercury = body;
                }
            }
            double theta_p_100, theta_p_0;
            theta_p_0 = atan2(mercury->r0(1), mercury->r0(0)); // Just zero

            double r_peri = 0.3075;
            double v_peri = 12.44;
            double current_r, current_v;
            double eps = 5e-8; // our precision should be at best ~5.e-9 (T=100, N=1e7). The result from using this epsilon was also pretty stable

            // Finding the last perihelion by iterating from the end
            for (int i = (N-1); i >= 0; i--) {
                current_r = norm(mercury->r.row(i),2);
                current_v = norm(mercury->v.row(i),2);
                if (abs(current_r - r_peri) < eps && abs(current_v - v_peri) < eps) {
                    theta_p_100 = atan2(mercury->r.row(i)(1), mercury->r.row(i)(0));
                    cout << "Found at iteration " << i << endl; // Important that I find it in the last orbit..
                    break;
                }
            }

            double precession_per_century = (theta_p_100 - theta_p_0); // radians
            cout << precession_per_century << " radians" << endl;
            cout << "Observed precession: " << precession_per_century*206300. << " arc seconds per century." << endl;

        }

    private:
        rowvec acceleration_gravity(rowvec r, rowvec r_other, double mass_solar) {
            rowvec R_vec = r - r_other;
            double R = norm(R_vec, 2);
            return -G*mass_solar*R_vec / pow(R, beta+1);
        }

        rowvec acceleration_mercury_GR(rowvec r, rowvec v) {
            double R = norm(r,2);
            double c, l;
            c = 63239.7263;    // AU/yr
            l = abs(r(0)*v(1) - r(1)*v(0));
            double GR = 1 + 3.*pow(l,2) / (pow(R,2)*pow(c,2));
            return -(G / pow(R,3)) * GR * r;
        }
};

int main() {
    // Compile as "c++ project3.cpp -larmadillo -o project3.exe", then run the executable

    int dim = 2;        // Spatial dimensions
    double T = 1.;      // yrs
    int N = 1000;
    double h = T / double(N);       // h = dt in yrs
    rowvec r0_earth = {1,0};
    rowvec v0_earth = {0,2*M_PI};   // [AU/yrs]. Circular
    rowvec origin = zeros<rowvec>(dim);
    double sun_mass = 2.E30;    // kg
    double earth_mass = 6.E24;  // kg

    ODESolver solver(dim, N, h);
    solver.add_body("Earth", earth_mass, r0_earth, v0_earth);
    solver.add_body("Sun", sun_mass, origin, origin);
    solver.set_inertial_sun(true);  // The acceleration of the Sun is now always zero. If the initial velocity is zero, then it is stationary through the whole simulation.

    // Let the code above stay un-commented

    double N_values[] = {1e2, 1e3, 1e4, 1e5, 1e6, 1e7};  // DON'T INCREASE ABOVE 1E7 (crashes my computer + memory error)
    double log_h, log_error;
    clock_t start, finish;
    double time_used;


    // // Euler errors!
    // ofstream outfile1;
    // outfile1.open("errors_euler.txt");
    // for (int n : N_values) {
    //     h = T / double(n);
    //     solver.set_resolution(n,h);
    //     start = clock();
    //     solver.euler();
    //     finish = clock();
    //     time_used = (double) (finish - start) / CLOCKS_PER_SEC;
    //     cout << time_used << endl;
    //     log_error = solver.check_stability_circular();
    //     log_h = log10(h);
    //     // Write errors to file. Format: "log10(h),log10(eps)\n"
    //     outfile1 << setprecision(8) << log_h << "," << log_error << endl;
    // }
    // outfile1.close();

    // // Velocity Verlet errors!
    // ofstream outfile2;
    // outfile2.open("errors_verlet.txt");
    // for (int n : N_values) {
    //     h = T / double(n);
    //     solver.set_resolution(n,h);
    //     start = clock();
    //     solver.velocity_verlet();
    //     finish = clock();
    //     time_used = (double) (finish - start) / CLOCKS_PER_SEC;
    //     cout << time_used << endl;
    //     log_error = solver.check_stability_circular();
    //     log_h = log10(h);
    //     // Write errors to file. Format: "log10(h),log10(eps)\n"
    //     outfile2 << setprecision(8) << log_h << "," << log_error << endl;
    // }
    // outfile2.close();



    // Optimal Verlet
    // N = 1e7;
    // h = T / double(N);
    // solver.set_resolution(N, h);
    // solver.velocity_verlet();
    // solver.write_to_file("Earth", "optimal_verlet.dat");


    // // Print errors for N=1e7
    // bool print = true;
    // cout << "Euler:\n";
    // solver.euler();
    // solver.check_stability_circular(print);
    // cout << endl;
    // cout << "Verlet:\n";
    // solver.velocity_verlet();
    // solver.check_stability_circular(print);


    // // Escape velocity testing
    // T = 1000.;
    // N = 1e7;
    // h = T / double(N);
    // solver.set_resolution(N,h);
    // rowvec v_escape_analytic = {0,2*sqrt(2)*M_PI};     // 8.8858 AU/yr
    // rowvec v_escape_attempt = {0,8.8};
    // bool bound = true;

    // while (bound) {
    //     v_escape_attempt(1) += .02;
    //     solver.initialize_body("Earth", r0_earth, v_escape_attempt);
    //     solver.velocity_verlet();
    //     bound = solver.check_if_enclosed();
    // }
    // solver.write_to_file("Earth", "escape.dat");
    // cout << "Analytical escape velocity: " << norm(v_escape_analytic,2) << "AU/yr." << endl;
    // cout << "Numerical escapce velocity: " << norm(v_escape_attempt,2) << "AU/yr." << endl;


    // // Beta testing
    // T = 5;
    // N = 1e7;
    // h = T / double(N);
    // solver.set_resolution(N,h);
    // rowvec v0_beta = {0,8.5};
    // solver.initialize_body("Earth", r0_earth, v0_earth);
    // double beta[] = {3.};//{2., 2.2, 2.4, 2.6, 2.8, 3.};
    // string filename;
    // for (double b : beta) {
    //     solver.setBeta(b);
    //     solver.velocity_verlet();
    //     filename = "beta" + to_string(b) + ".dat";
    //     solver.write_to_file("Earth", filename);
    // }


    // Adding Jupiter!
    // T = 13.;
    // N = 1e6;
    // h = T / double(N);
    // double jupiter_mass = 1.9E27;
    // solver.set_resolution(N, h);
    // rowvec r0_jup = {5.2, 0};
    // rowvec v0_jup = {0,2*M_PI / sqrt(5.2)};
    // solver.add_body("Jupiter", jupiter_mass, r0_jup, v0_jup);
    // solver.initialize_body("Earth", r0_earth, v0_earth);
    // solver.set_body_mass("Jupiter", jupiter_mass);
    // solver.velocity_verlet();
    // solver.write_to_file();
    // solver.check_stability_circular(true);



    // // Allow the sun to move!
    // solver.set_inertial_sun(false);
    // // Finding new start positions in cm system.
    // double m_tot = jupiter_mass + sun_mass + earth_mass;
    // rowvec R_cm = pow(m_tot,-1) * (earth_mass*r0_earth + jupiter_mass*r0_jup + sun_mass*origin);
    // rowvec r0_sun = -R_cm;
    // r0_earth -= R_cm;
    // r0_jup -= R_cm;
    // rowvec v0_sun = pow(sun_mass,-1)*( -earth_mass*v0_earth - jupiter_mass*v0_jup); // Total momentum should be zero, so cm is fixed
    //
    // solver.initialize_body("Earth", r0_earth, v0_earth);
    // solver.initialize_body("Jupiter", r0_jup, v0_jup);
    // solver.initialize_body("Sun", r0_sun, v0_sun);
    //
    // solver.velocity_verlet();
    // // solver.write_to_file();
    // solver.check_stability_circular(true);



    // // Finding the orbits of all planet
    // T = 1.;
    // N = 1e6;
    // h = T / double(N);
    // dim = 3;
    // ODESolver solver_3D(dim, N, h);
    // solver_3D.set_inertial_sun(false); // Allow star to move
    // void add_all_planets(ODESolver*);
    // add_all_planets(&solver_3D); // + Sun
    // solver_3D.velocity_verlet();
    // solver_3D.write_to_file();



    // Studying the precession of Mercury
    dim = 2;
    T = 100.;
    N = 1e7;
    h = T / double(N);
    ODESolver merc_orbit(dim, N, h);
    merc_orbit.set_inertial_sun(true);  // "orbit of Mercury around the Sun"
    rowvec r0_merc = {0.3075,0};        // Perihelion
    rowvec v0_merc = {0,12.44};         // Corresponding velocity at perihelion
    double mass_merc = 3.3E23;
    merc_orbit.add_body("Mercury", mass_merc, r0_merc, v0_merc);
    merc_orbit.add_body("Sun", sun_mass, origin, origin);
    merc_orbit.velocity_verlet();
    cout << "Without GR term: \n";
    merc_orbit.calculate_precession(); // Without GR term
    cout << endl;

    merc_orbit.set_mercury_GR(true);
    merc_orbit.velocity_verlet();
    cout << "With GR term: \n";
    merc_orbit.calculate_precession();
    // merc_orbit.write_to_file("Mercury", "merc_GR.dat");


    return 0;
}

void add_all_planets(ODESolver* slvr) {
    // Mars: 2458410.500000000 = A.D. 2018-Oct-19 00:00:00.0000 TDB (NASA)
    double mars_mass = 6.6E23;
    rowvec r0_mars = {1.386064640638050E+00,-7.195522876836125E-02,-3.574720875419098E-02};
    rowvec v0_mars = {1.325375632613797E-03,1.516910673406354E-02, 2.852845787253853E-04};  // AU / day
    v0_mars *=  365.; // AU / yr
    slvr->add_body("Mars", mars_mass, r0_mars, v0_mars);

    double sun_mass = 2.E30;
    rowvec r0_sun = {-1.854398580107027E-04,7.261765071449882E-03,-7.173229132111148E-05};
    rowvec v0_sun = {-7.598223688385497E-06,2.546589792171653E-06,1.899439973348268E-07};
    v0_sun *= 365.;
    slvr->add_body("Sun", sun_mass, r0_sun, v0_sun);

    double mercury_mass = 3.3E23;
    rowvec r0_merc = {-4.466974942983433E-02,-4.551297184815911E-01,-3.377443034243644E-02};
    rowvec v0_merc = {2.235388667387174E-02,-1.255928387575858E-03,-2.154047309352793E-03};
    v0_merc *= 365.;
    slvr->add_body("Mercury", mercury_mass, r0_merc, v0_merc);

    double venus_mass = 4.9E24;
    rowvec r0_venus = {6.771183423096696E-01,2.635570892119448E-01,-3.564015770708658E-02};
    rowvec v0_venus = {-7.233924240439758E-03,1.883053303306194E-02,6.755080544414482E-04};
    v0_venus *= 365;
    slvr->add_body("Venus", venus_mass, r0_venus, v0_venus);

    double earth_mass =  6.E24;
    rowvec r0_earth = {9.004267194046488E-01,4.329702891250327E-01,-9.309259935529284E-05};
    rowvec v0_earth = {-7.644048786784979E-03,1.548720517812966E-02,9.799447004349826E-08};
    v0_earth *= 365.;
    slvr->add_body("Earth", earth_mass, r0_earth, v0_earth);

    double jupiter_mass;
    rowvec r0_jupiter = {-2.628114780107304E+00,-4.675856649493857E+00,7.818020939743552E-02};
    rowvec v0_jupiter = {6.487854045762937E-03,-3.337752963125740E-03,-1.312915365947336E-04};
    v0_jupiter *= 365;
    slvr->add_body("Jupiter", jupiter_mass, r0_jupiter, v0_jupiter);

    double saturn_mass = 5.5E26;
    rowvec r0_saturn = {1.575176167525727E+00,-9.930952889722587E+00,1.099718612784312E-01};
    rowvec v0_saturn = {5.202593393091781E-03,8.555506892745189E-04,-2.216362261569332E-04};
    v0_saturn *= 365;
    slvr->add_body("Saturn", saturn_mass, r0_saturn, v0_saturn);

    double uranus_mass = 8.8E25;
    rowvec r0_uranus = {1.716586975644875E+01,1.001373868288497E+01,-1.851949198554858E-01};
    rowvec v0_uranus = {-2.010719471096281E-03,3.213990583409865E-03,3.805302445255308E-05};
    v0_uranus *= 365;
    slvr->add_body("Uranus", uranus_mass, r0_uranus, v0_uranus);

    double neptune_mass =  1.03E26;
    rowvec r0_neptune = {2.892424054587586E+01,-7.707281387885470E+00,-5.078715497631403E-01};
    rowvec v0_neptune = {7.870585771863339E-04,3.052076356482153E-03,-8.060978991971341E-05};
    v0_neptune *= 365;
    slvr->add_body("Neptune", neptune_mass, r0_neptune, v0_neptune);
}
