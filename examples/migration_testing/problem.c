/**
 * Obliquity Sculpting of Kepler Multis (Millholland & Laughlin 2019)
 *
 * In particular, this simulates the obliquity sculpting of the Kepler multis due to convergent migration over 4 Myr.
 * Two planets are initialized just wide of the 3:2 MMR (0.173 and 0.233 AU) are migrated inward. After 2 Myr, migration is turned off. After a period of chaotic obliquity evolution, the inner planet is excited to high obliquity and maintained indefinitely
 * Result is based on Figure 3 in Millholland & Laughlin (2019), and reproduces Figure 3 in Lu et. al (2023)
 * For a more in-depth description of the various parameters that can be set in this simulation, please see the ipython examples for consistent tides & spin (any notebook with the prefix TidesSpin)
 * and migration (Migration.ipynb)
 *
 * integration time is artificially shortened to run quickly. Full run in paper with tmax = 4e6 orbits takes ~ 2 days
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "tides_spin.c"

void heartbeat(struct reb_simulation* sim);
//double factor = 5.;
double tmax = 4e6*2*M_PI; // set short to run quickly. Set to 4e6 * 2 * M_PI in paper

int Ntest = 1;
double inv_alignment_ts = 1./(1e5 * M_PI * 2.);

char title1[100] = "output_orbits_1018.txt";
char title2[100] = "output_spins_1018.txt";
char title3[100] = "output_tp_1010_nd_";

const int p = 2;

int main(int argc, char* argv[]){

    int index = 0;
    if (argc == 2){
      strcat(title3, argv[1]);
      index = atoi(argv[1]);
    }

    struct reb_simulation* sim = reb_create_simulation();
    // Exact parameters from Millholland & Laughlin (2019)
    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    //reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object
    const double pm = 5. * 3.0e-6;
    const double pr = 2.5 * 4.26e-5;

    //reb_add_fmt(sim, "m a e r inc Omega pomega M", pm, 0.17308688, 0.01, pr, 0.5 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 1
    //reb_add_fmt(sim, "m a e r inc Omega pomega M", pm, 0.23290608, 0.01, pr, -0.431 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 2

    struct reb_particle star = {};
    star.m = solar_mass;
    star.r = solar_rad;
    star.x = -1.137884e-06;
    star.y = -5.255686e-06;
    star.z = -2.030614e-08;
    star.vx = 5.690636e-05;
    star.vy = -1.872905e-05;
    star.vz = -1.548650e-08;

    struct reb_particle p1 = {};
    p1.m = pm;
    p1.r = pr;
    p1.x = 1.203485e-01;
    p1.y = 1.239405e-01;
    p1.z = 1.100458e-03;
    p1.vx = -1.743498e+00;
    p1.vy = 1.660694e+00;
    p1.vz = -1.451128e-02;

    struct reb_particle p2 = {};
    p2.m = pm;
    p2.r = pr;
    p2.x = -4.448906e-02;
    p2.y = 2.264377e-01;
    p2.z = 2.532830e-04;
    p2.vx = -2.050259e+00;
    p2.vy = -4.120908e-01;
    p2.vz = 1.554372e-02;

    reb_add(sim, star);
    reb_add(sim, p1);
    reb_add(sim, p2);
    // add test particles
    sim->N_active = 3;
    sim->integrator = REB_INTEGRATOR_IAS15;
    //struct reb_orbit op1 = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    // sim->dt = 1e-3;//op1.P / 10.56789;
    sim->heartbeat = heartbeat;

    // Add REBOUNDx Additional effects
    // First Spin
    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    //struct rebx_force* damping = rebx_load_force(rebx, "laplace_damping");
    rebx_add_force(rebx, effect);
    //rebx_add_force(rebx, damping);
    // Exact parameters from Millholland & Laughlin (2019)
    // Sun
    const double solar_spin_period = 20 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_Q = 1000000.;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", 0.07);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * solar_mass * solar_rad * solar_rad);

    struct reb_vec3d solar_spin_vec = {};
    solar_spin_vec.x = 0.000000e+00;
    solar_spin_vec.y = 1.633457e-06;
    solar_spin_vec.z = 1.825000e+01;

    double orbn = 1.388696e+01;
    double orb2n = 8.895784e+00;

    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", solar_spin_vec);
    //rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z = solar_spin}); // Omega_x = Omega_y = 0 by default

    // We assume tau = 1/(2*n*Q) with n the mean motion, even though the spin is not synchronized with the orbit (see Lu et al. (2023))
    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", 1./(2.*orbn*solar_Q));

    // P1
    const double spin_period_1 = 5. * 2. * M_PI / 365.; // 5 days in REBOUND time units
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * pm * pr * pr);

    struct reb_vec3d p1_spin_vec = {};
    p1_spin_vec.x = 0.000000e+00;
    p1_spin_vec.y = -1.910907e+00;
    p1_spin_vec.z = 7.297498e+01;

    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", p1_spin_vec);
    //rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", (struct reb_vec3d){.y=spin_1 * -0.0261769, .z=spin_1 * 0.99965732});

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2.*orbn*planet_Q));

    // P2
    double spin_period_2 = 3. * 2. * M_PI / 365.; // 3 days in REBOUND time units
    double spin_2 = (2. * M_PI) / spin_period_2;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "I", 0.25 * pm * pr * pr);

    struct reb_vec3d p2_spin_vec = {};
    p2_spin_vec.x = 0.000000e+00;
    p2_spin_vec.y = 3.038466e+00;
    p2_spin_vec.z = 1.216287e+02;
    rebx_set_param_vec3d(rebx, &sim->particles[2].ap, "Omega", p2_spin_vec);
    //rebx_set_param_vec3d(rebx, &sim->particles[2].ap, "Omega", (struct reb_vec3d){.y=spin_2 * 0.0249736, .z=spin_2 * 0.99968811});

    //struct reb_orbit orb2 = reb_tools_particle_to_orbit(sim->G, sim->particles[2], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau", 1./(2.*orb2n*planet_Q));

    struct reb_particle planet = sim->particles[p];
    struct reb_vec3d* Omega_p = rebx_get_param(rebx, planet.ap, "Omega");

    // Rotation from inv frame to planet frame
    struct reb_orbit orbp = reb_tools_particle_to_orbit(sim->G, planet, sim->particles[0]);
    struct reb_vec3d lonp = reb_vec3d_cross((struct reb_vec3d){.z =1}, orbp.hvec);  // Line of nodes is the new x-axis
    struct reb_rotation rotp = reb_rotation_init_to_new_axes(orbp.hvec, lonp);

    struct reb_vec3d svp = reb_vec3d_rotate(*Omega_p, rotp);

    double magp;
    double thetap;
    double phip;
    reb_tools_xyz_to_spherical(svp, &magp, &thetap, &phip);

    double maginv;
    double thetainv;
    double phiinv;
    reb_tools_xyz_to_spherical(*Omega_p, &maginv, &thetainv, &phiinv);

    double mag_orb;
    double theta_orb;
    double phi_orb;
    reb_tools_xyz_to_spherical(orbp.hvec, &mag_orb, &theta_orb, &phi_orb);

    // Test particles
    const double j2 = (0.4 / 3.) * (magp * magp * pr * pr * pr) / (sim->G * pm);
    const double lr = pow(2. * j2 * pr * pr * orbp.a * orbp.a * orbp.a * pow((1 - orbp.e),(3./2.)) * pm / solar_mass, (1./5.));
    //for (unsigned int i = 0; i < Ntest; i++){

    double d = (1. + 0.2 * index) * pr;
    double le_theta = thetap - 0.5 * atan(sin(2.*thetap)/(cos(2.*thetap) + (lr/d)*(lr/d)*(lr/d)*(lr/d)*(lr/d)));
    //printf("%f\n", le_theta * 180./M_PI);
    // Initialize particles on the Laplace surface!
    // Initialize in the planet frame
    //reb_add_fmt(sim, "primary a inc Omega", planet, d, theta_orb + le_theta, 90 * M_PI/180. + phiinv);
/*
    struct reb_orbit ob = reb_tools_particle_to_orbit(sim->G, sim->particles[3], sim->particles[p]);
    double magb;
    double thetab;
    double phib;
    reb_tools_xyz_to_spherical(ob.hvec, &magb, &thetab, &phib);
*/
      //printf("TP Before rot: %f %f\n", thetab*180./M_PI, phib*180./M_PI);
      // reb_particle_irotate(&sim->particles[3], rotp);
    //}
    //exit(1);

    // And migration
    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    // Set migration parameters
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -5e6 * 2 * M_PI);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", (-5e6 * 2 * M_PI) / 1.1);

    reb_move_to_com(sim);

    // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
    // nodes with the invariable plane (along z cross newz)
    //struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    //struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    //struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    //rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
    rebx_spin_initialize_ode(rebx, effect);

    // Run simulation
    //system("rm -v output_orbits.txt"); // remove previous output files
    //system("rm -v output_spins.txt");
    //system("rm -v output_tp.txt");

    FILE* of_orb = fopen(title1, "a");
    fprintf(of_orb,"t,a1,Om1,a2,Om2\n");
    fclose(of_orb);

    FILE* of_spins = fopen(title2, "a");
    fprintf(of_spins,"t,mag1,theta1,mag2,theta2\n");
    fclose(of_spins);
/*
    FILE* of = fopen(title3, "w");
    fprintf(of, "t");
    fprintf(of, ",nx1,ny1,nz1,a1,theta1,phi1");
    fprintf(of,",lr\n");
    fclose(of);
*/
    //system("rm -v output_test_1.txt");

    reb_integrate(sim, tmax/2);
    //tmax *= factor;
    //reb_integrate(sim,3tmax);
    // reb_simulationarchive_snapshot(sim, "archive.bin");

    //printf("Migration Switching Off\n");
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", INFINITY);

    reb_integrate(sim, tmax);

    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
  //if(reb_output_check(sim, tmax/100000)){        // outputs every 100 REBOUND years
  if(reb_output_check(sim, 10. * 2. * M_PI)){
    struct rebx_extras* const rebx = sim->extras;
    // FILE* of_test = fopen(title3, "a");
    //if (of_orb == NULL || of_spins == NULL){
    //    reb_error(sim, "Can not open file.");
    //    return;
    //}
    struct reb_particle* sun = &sim->particles[0];

    FILE* of_orb = fopen(title1, "a");
    FILE* of_spins = fopen(title2, "a");
    struct reb_particle* p1 = &sim->particles[1];
    struct reb_particle* p2 = &sim->particles[2];

    struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
    double a1 = o1.a;
    double e1 = o1.e;
    double i1 = o1.inc;
    double Om1 = o1.Omega;
    double pom1 = o1.pomega;
    struct reb_vec3d norm1 = o1.hvec;

    struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *p2, *sun);
    double a2 = o2.a;
    double e2 = o2.e;
    double i2 = o2.inc;
    double Om2 = o2.Omega;
    double pom2 = o2.pomega;
    struct reb_vec3d norm2 = o2.hvec;

    // Spin vectors - all initially in invariant plane
    struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

    // Interpret both planet spin vectors in the rotating planet frame in spherical coordinates
    struct reb_vec3d* Omega_p1 = rebx_get_param(rebx, p1->ap, "Omega");

    struct reb_vec3d lon1 = reb_vec3d_cross((struct reb_vec3d){.z =1}, norm1);  // Line of nodes is the new x-axis
    struct reb_rotation rot1 = reb_rotation_init_to_new_axes(norm1, lon1);      // Arguments to this function are the new z and x axes
    struct reb_vec3d sv1 = reb_vec3d_rotate(*Omega_p1, rot1);

    double mag1;
    double theta1;
    double phi1;
    reb_tools_xyz_to_spherical(sv1, &mag1, &theta1, &phi1);

    struct reb_vec3d* Omega_p2 = rebx_get_param(rebx, p2->ap, "Omega");

    struct reb_vec3d lon2 = reb_vec3d_cross((struct reb_vec3d){.z =1}, norm2); // Line of nodes is the new x-axis
    struct reb_rotation rot2 = reb_rotation_init_to_new_axes(norm2, lon2); // Arguments to this function are the new z and x axes
    struct reb_vec3d sv2 = reb_vec3d_rotate(*Omega_p2, rot2);

    double mag2;
    double theta2;
    double phi2;
    reb_tools_xyz_to_spherical(sv2, &mag2, &theta2, &phi2);



    //fprintf(of_orb, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, a1, e1, i1, pom1, Om1, norm1.x, norm1.y, norm1.z, a2, e2, i2, pom2, Om2, norm2.x, norm2.y, norm2.z);  // prints the spins and orbits of all bodies
    fprintf(of_orb, "%e,%e,%e,%e,%e\n", sim->t, a1, Om1, a2, Om2);  // prints the spins and orbits of all bodies
    fclose(of_orb);
    //fprintf(of_spins, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, Omega_sun->x, Omega_sun->y, Omega_sun->z, mag1, theta1, phi1, mag2, theta2, phi2, Omega_p1->x, Omega_p1->y, Omega_p1->z, Omega_p2->x, Omega_p2->y, Omega_p2->z);
    fprintf(of_spins, "%e,%e,%e,%e,%e\n", sim->t, mag1, theta1, mag2, theta2);
    fclose(of_spins);
/*
    struct reb_particle* planet = &sim->particles[p];
    // Orbit information

    struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *planet, *sun);
    double ap = o.a;
    double ep = o.e;
    double ip = o.inc;
    double Omp = o.Omega;
    double pomp = o.pomega;
    struct reb_vec3d normp = o.hvec;

    struct reb_vec3d lon = reb_vec3d_cross((struct reb_vec3d){.z =1}, normp);  // Line of nodes is the new x-axis
    struct reb_rotation rot = reb_rotation_init_to_new_axes(normp, lon);      // Arguments to this function are the new z and x axes

    struct reb_vec3d* Omega_p = rebx_get_param(rebx, planet->ap, "Omega");
    struct reb_vec3d svp = reb_vec3d_rotate(*Omega_p, rot);

    double magp;
    double thetap;
    double phip;
    reb_tools_xyz_to_spherical(svp, &magp, &thetap, &phip);

    struct reb_particle* test = &sim->particles[3];
    struct reb_orbit ot = reb_tools_particle_to_orbit(sim->G, *test, *planet);
    double at = ot.a;
    double et = ot.e;
    double it = ot.inc;
    double Omt = ot.Omega;
    double pomt = ot.pomega;
    struct reb_vec3d normt = ot.hvec;

    double mag;
    double theta;
    double phi;

    reb_vec3d_irotate(&normt, rot); // rotates into planet frame
    reb_tools_xyz_to_spherical(normt, &mag, &theta, &phi);


    const double* k2 = rebx_get_param(sim->extras, sim->particles[p].ap, "k2");
    double rp = sim->particles[p].r;
    double mp = sim->particles[p].m;
    double ms = sim->particles[0].m;
    const double j2 = ((*k2) / 3.) * (magp * magp * rp * rp * rp) / (sim->G * mp);
    const double lr = pow(2. * j2 * rp * rp * ap * ap * ap * pow((1 - ep * ep),(3./2.)) * mp / ms, (1./5.));
    //printf("%f\n",lr);
    fprintf(of_test, "%f,%f,%f,%f,%f,%f,%f,%f\n", sim->t, normt.x, normt.y, normt.z, at, theta, phi,lr);
    fclose(of_test);
*/
  }

  if(reb_output_check(sim, 100.*M_PI)){        // outputs to the screen
      reb_output_timing(sim, tmax);
  }
}
