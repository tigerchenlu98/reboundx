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
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "tides_spin.c"

void heartbeat(struct reb_simulation* sim);
double factor = 5.;
double tmax = 1e6; // set short to run quickly. Set to 4e6 * 2 * M_PI in paper

int Ntest = 10;
double inv_alignment_ts = 1./(1e6 * M_PI * 2.);

void derivatives(struct reb_ode* const ode, double* const yDot, const double* const y, const double t){
    // From Zanazzi & Lai 2017
    // hard-coded relevant physical parameters
    const double ms = ode->r->particles[0].m;
    const double mp = ode->r->particles[1].m;
    const double rp = ode->r->particles[1].r;
    struct reb_orbit op = reb_tools_particle_to_orbit(ode->r->G, ode->r->particles[1], ode->r->particles[0]); // planet orbit

    struct reb_vec3d lp = op.hvec;
    struct reb_vec3d lp_hat = reb_vec3d_normalize(lp);

    const struct reb_vec3d* const sp_ptr = rebx_get_param(ode->r->extras, ode->r->particles[1].ap, "Omega");
    struct reb_vec3d sp;
    if (sp_ptr != NULL){
      sp = *sp_ptr;
    }
    struct reb_vec3d sp_hat = reb_vec3d_normalize(sp);

    // This stuff for alignment torque
    struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, op.hvec);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(op.hvec, line_of_nodes); // Invariant plane to planet
    if (isnan(rot.r)) {
      rot = reb_rotation_identity();
    }
    struct reb_vec3d sp_rot = reb_vec3d_rotate(sp, rot);

    double omega;
    double theta_p;
    double phi_p;
    reb_tools_xyz_to_spherical(sp_rot, &omega, &theta_p, &phi_p); // theta and phi in the planet frame
    struct reb_rotation invrot = reb_rotation_inverse(rot); // planet frame to invariant plane

    // Planet quadrupole
    const double* k2 = rebx_get_param(ode->r->extras, ode->r->particles[1].ap, "k2");
    const double j2 = ((*k2) / 3.) * (omega * omega * rp * rp * rp) / (ode->r->G * mp);
    const double lr = pow(j2 * rp * rp * op.a * op.a * op.a * pow((1 - op.e * op.e),(3./2.)) * mp / ms, (1./5.));

    for (unsigned int i = 0; i < Ntest; i++){
      // Multiple particles
      struct reb_vec3d l_hat = {};
      l_hat.x = y[i*4+0];
      l_hat.y = y[i*4+1];
      l_hat.z = y[i*4+2];

      // Planet-TP distance
      double d = y[i*4+3];

      // star
      const double star_prefactor = ((3. * ode->r->G * ms * d * d) / (4 * op.a * op.a * op.a)) * reb_vec3d_dot(l_hat, lp_hat);;
      struct reb_vec3d l_cross_lp = reb_vec3d_cross(l_hat, lp_hat);
      struct reb_vec3d tstar = reb_vec3d_mul(l_cross_lp, star_prefactor);

      // planet
      const double planet_prefactor = ((3. * ode->r->G * mp * rp * rp * j2) / (2. * d * d * d)) * reb_vec3d_dot(l_hat, sp_hat);
      struct reb_vec3d l_cross_sp = reb_vec3d_cross(l_hat, sp_hat);
      struct reb_vec3d tplanet = reb_vec3d_mul(l_cross_sp, planet_prefactor);

      // simple alignment torque
      // Align towards Laplace Equilibrium in planet frame
      const double le_theta = theta_p - 0.5 * atan2(sin(2.*theta_p),cos(2.*theta_p) + (lr/d)*(lr/d)*(lr/d)*(lr/d)*(lr/d)); // In the planet frame
      struct reb_vec3d beta_vec = {};
      beta_vec.z = 1.;

      //struct reb_vec3d beta_vec = reb_tools_spherical_to_xyz(1., le_theta, phi_p);// In the planet frame
      reb_vec3d_irotate(&beta_vec, invrot); // rotates into xyz frame.

      struct reb_vec3d beta_cross_l = reb_vec3d_cross(beta_vec, l_hat);
      struct reb_vec3d l_cross_t1 = reb_vec3d_cross(l_hat, beta_cross_l);
      struct reb_vec3d talign = reb_vec3d_mul(l_cross_t1, 1.);//inv_alignment_ts);

      // DiffEq
      const double inv_prefactor = 1. / (d * d * sqrt(ode->r->G * mp / (d * d * d)));
      yDot[i*4+0] = inv_prefactor * (tstar.x + tplanet.x) + talign.x;
      yDot[i*4+1] = inv_prefactor * (tstar.y + tplanet.y) + talign.y;
      yDot[i*4+2] = inv_prefactor * (tstar.z + tplanet.z) + talign.z;
      yDot[i*4+3] = 0;
    }
}

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Exact parameters from Millholland & Laughlin (2019)
    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    const double pm = 5. * 3.0e-6;
    const double pr = 2.5 * 4.26e-5;
    //const double pa1_init = 0.5;
    //const double pa2_init = pa1_init * pow(1.55, 2./3.);
    //reb_add_fmt(sim, "m r a inc", pm, pr, pa1_init, 0.3 * (M_PI/180.)); // Planet 1
    //reb_add_fmt(sim, "m r a inc", pm, pr, pa2_init, -0.3 * (M_PI/180.)); // Planet 1

    reb_add_fmt(sim, "m a e r inc Omega pomega M", pm, 0.17308688, 0.01, pr, 0.5 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 1
    reb_add_fmt(sim, "m a e r inc Omega pomega M", pm, 0.23290608, 0.01, pr, -0.431 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 2

    // add test particles
    sim->N_active = 3;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    //struct reb_orbit op1 = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    sim->dt = 1e-3;//op1.P / 10.56789;
    sim->heartbeat = heartbeat;

    struct reb_particle planet = sim->particles[1];

    struct reb_ode* ho = reb_create_ode(sim,Ntest*4+1);   // Add an ODE with 2 dimensions
    ho->derivatives = derivatives;              // Right hand side of the ODE

    // Test particle
    for (unsigned int i = 0; i < Ntest; i++){
      double ta = (1.5 + (double)i * 0.5) * pr;
      reb_add_fmt(sim, "primary a", planet, ta);

      struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[2], sim->particles[1]);
      struct reb_vec3d l = o.hvec;
      struct reb_vec3d l_hat = reb_vec3d_normalize(l);

      //if (i == 0){
      //  sim->dt = o.P / 10.6789;
      //}

      ho->y[i*4+0] = l_hat.x;                               // Initial conditions
      ho->y[i*4+1] = l_hat.y;
      ho->y[i*4+2] = l_hat.z;
      ho->y[i*4+3] = ta;

      // remove particle
      reb_remove(sim, 3, 1);
    }

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
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z = solar_spin}); // Omega_x = Omega_y = 0 by default

    // We assume tau = 1/(2*n*Q) with n the mean motion, even though the spin is not synchronized with the orbit (see Lu et al. (2023))
    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", 1./(2.*orb.n*solar_Q));

    // P1
    const double spin_period_1 = 5. * 2. * M_PI / 365.; // 5 days in REBOUND time units
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * pm * pr * pr);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", (struct reb_vec3d){.y=spin_1 * -0.0261769, .z=spin_1 * 0.99965732});

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2.*orb.n*planet_Q));

    // P2
    double spin_period_2 = 3. * 2. * M_PI / 365.; // 3 days in REBOUND time units
    double spin_2 = (2. * M_PI) / spin_period_2;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "I", 0.25 * pm * pr * pr);
    rebx_set_param_vec3d(rebx, &sim->particles[2].ap, "Omega", (struct reb_vec3d){.y=spin_2 * 0.0249736, .z=spin_2 * 0.99968811});

    struct reb_orbit orb2 = reb_tools_particle_to_orbit(sim->G, sim->particles[2], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau", 1./(2.*orb2.n*planet_Q));

    // And migration
    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    // Set migration parameters
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -5e6 * 2 * M_PI);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", (-5e6 * 2 * M_PI) / 1.1);

    reb_move_to_com(sim);

    // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
    // nodes with the invariable plane (along z cross newz)
    struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
    rebx_spin_initialize_ode(rebx, effect);

    // Run simulation
    //system("rm -v output_orbits.txt"); // remove previous output files
    //system("rm -v output_spins.txt");
    //system("rm -v output_tp.txt");

    FILE* of_orb = fopen("output_orbits_ml19_6.txt", "a");
    fprintf(of_orb,"t,a1,Om1,a2,Om2\n");
    fclose(of_orb);

    FILE* of_spins = fopen("output_spins_ml19_6.txt", "a");
    fprintf(of_spins,"t,mag1,theta1,mag2,theta2\n");
    fclose(of_spins);

    FILE* of = fopen("output_tp_ml19_bvec_6.txt", "w");
    fprintf(of, "t");
    for (int i = 1; i < Ntest+1; i++){
      fprintf(of, ",nx%d,ny%d,nz%d,a%d,theta%d,phi%d",i,i,i,i,i,i);
    }
    fprintf(of,",lr\n");
    fclose(of);

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
  if(reb_output_check(sim, tmax/100000)){        // outputs every 100 REBOUND years
    struct rebx_extras* const rebx = sim->extras;
    FILE* of_orb = fopen("output_orbits_ml19_6.txt", "a");
    FILE* of_spins = fopen("output_spins_ml19_6.txt", "a");
    FILE* of_test = fopen("output_tp_ml19_bvec_6.txt", "a");
    //if (of_orb == NULL || of_spins == NULL){
    //    reb_error(sim, "Can not open file.");
    //    return;
    //}

    struct reb_particle* sun = &sim->particles[0];
    struct reb_particle* p1 = &sim->particles[1];
    struct reb_particle* p2 = &sim->particles[2];

    // Orbit information

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

/*
    struct reb_orbit ot = reb_tools_particle_to_orbit(sim->G, *test, *p1);
    double at = ot.a;
    double et = ot.e;
    double it = ot.inc;
    double Omt = ot.Omega;
    double pomt = ot.pomega;
    struct reb_vec3d normt = ot.hvec;
*/
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
    //fprintf(of_test, "%f,%f,%f,%f,%f,%f,%f\n", sim->t, at, et, it, normt.x, normt.y, normt.z);
    //fclose(of_test);

    struct reb_ode** odes = sim->odes;
    struct reb_ode* lhat_ode = odes[0];

    fprintf(of_test, "%f",sim->t);

    for (unsigned int i = 0; i < Ntest; i++){
      struct reb_vec3d lhat = {};
      lhat.x = lhat_ode->y[i*4+0];
      lhat.y = lhat_ode->y[i*4+1];
      lhat.z = lhat_ode->y[i*4+2];

      double mag;
      double theta;
      double phi;

      reb_vec3d_irotate(&lhat, rot1); // rotates into planet frame
      reb_tools_xyz_to_spherical(lhat, &mag, &theta, &phi);

      fprintf(of_test, ",%e,%e,%e,%e,%f,%f", lhat.x, lhat.y, lhat.z, lhat_ode->y[i*4+3],theta,phi);
    }

    const double* k2 = rebx_get_param(sim->extras, sim->particles[1].ap, "k2");
    double rp = sim->particles[1].r;
    double mp = sim->particles[1].m;
    double ms = sim->particles[0].m;
    const double j2 = ((*k2) / 3.) * (mag1 * mag1 * rp * rp * rp) / (sim->G * mp);
    const double lr = pow(j2 * rp * rp * o1.a * o1.a * o1.a * pow((1 - o1.e * o1.e),(3./2.)) * mp / ms, (1./5.));

    fprintf(of_test, ",%e\n", lr);
    fclose(of_test);
  }

  if(reb_output_check(sim, 100.*M_PI)){        // outputs to the screen
      reb_output_timing(sim, 6.*tmax);
  }
}
