/**
 * Self-Consistent Spin-Tidal and Dynamical equations of motion (Eggleton et. al 1998)
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
double tmax = 500 * 2 * M_PI; // set short to run quickly. Set to 4e6 * 2 * M_PI in paper

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Exact parameters from Millholland & Laughlin (2019)
    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    const double p1_mass = 5. * 3.0e-6; // in Earth masses * 1 Earth Mass / 1 Solar Mass
    const double p1_rad = 2.5 * 4.26e-5; // in Earth rad * 1 Earth rad / 1 AU
    reb_add_fmt(sim, "m a e r inc Omega pomega M", p1_mass, 0.17308688, 0.01, p1_rad, 0.5 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 1

    const double p2_mass = 5. * 3.0e-6;
    const double p2_rad = 2.5 * 4.26e-5;
    reb_add_fmt(sim, "m a e r inc Omega pomega M", p2_mass, 0.23290608, 0.01, p2_rad, -0.431 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 2
    sim->N_active = 3;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = 1e-3;
    sim->heartbeat = heartbeat;

    // Add REBOUNDx Additional effects
    // First Spin
    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);
    // Sun
    const double solar_spin_period = 20 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_Q = 1000000.;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", 0.07);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "moi", 0.07 * solar_mass * solar_rad * solar_rad);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sx", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sy", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sz", solar_spin * 1.0);
    double solar_sigma = rebx_tides_calc_sigma_from_Q(rebx, &sim->particles[0], &sim->particles[1], solar_Q);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sigma", solar_sigma);

    // P1
    const double spin_period_1 = 5. * 2. * M_PI / 365.; // 5 days in REBOUND time units
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "moi", 0.25 * p1_mass * p1_rad * p1_rad);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sx", spin_1 * 0.0);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sy", spin_1 * -0.0261769);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sz", spin_1 * 0.99965732);
    double planet_sigma_1 = rebx_tides_calc_sigma_from_Q(rebx, &sim->particles[1], &sim->particles[0], planet_Q);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sigma", planet_sigma_1);

    // P2
    double spin_period_2 = 3. * 2. * M_PI / 365.; // 3 days in REBOUND time units
    double spin_2 = (2. * M_PI) / spin_period_2;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "moi", 0.25 * p2_mass * p2_rad * p2_rad);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "sx", spin_2 * 0.0);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "sy", spin_2 * 0.0249736);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "sz", spin_2 * 0.99968811);
    double planet_sigma_2 = rebx_tides_calc_sigma_from_Q(rebx, &sim->particles[2], &sim->particles[0], planet_Q);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "sigma", planet_sigma_2);

    // And migration
    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    // Set migration parameters
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -5e6 * 2 * M_PI);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", (-5e6 * 2 * M_PI) / 1.1);

    // printf("Are we even initializing");
    reb_move_to_com(sim);
    rebx_align_simulation(rebx);
    rebx_spin_initialize_ode(rebx, effect);

    // Run simulation
    system("rm -v output.txt"); // remove previous output file
    reb_integrate(sim, tmax/2);

    printf("Migration Switching Off\n");
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", INFINITY);
    
    reb_integrate(sim, tmax);

    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
  if(reb_output_check(sim, tmax/10)){        // outputs every 100 REBOUND years
    struct rebx_extras* const rebx = sim->extras;
    FILE* of = fopen("output.txt", "a");
    if (of==NULL){
        reb_error(sim, "Can not open file.");
        return;
    }

    struct reb_particle* sun = &sim->particles[0];
    struct reb_particle* p1 = &sim->particles[1];
    struct reb_particle* p2 = &sim->particles[2];

    double* sx_sun = rebx_get_param(rebx, sun->ap, "sx");
    double* sy_sun = rebx_get_param(rebx, sun->ap, "sy");
    double* sz_sun = rebx_get_param(rebx, sun->ap, "sz");
    double mag_sun = sqrt(*sx_sun * *sx_sun + *sy_sun * *sy_sun + *sz_sun * *sz_sun);

    double* sx_p1 = rebx_get_param(rebx, p1->ap, "sx");
    double* sy_p1 = rebx_get_param(rebx, p1->ap, "sy");
    double* sz_p1 = rebx_get_param(rebx, p1->ap, "sz");
    double mag_p1 = sqrt(*sx_p1 * *sx_p1 + *sy_p1 * *sy_p1 + *sz_p1 * *sz_p1);

    double* sx_p2 = rebx_get_param(rebx, p2->ap, "sx");
    double* sy_p2 = rebx_get_param(rebx, p2->ap, "sy");
    double* sz_p2 = rebx_get_param(rebx, p2->ap, "sz");
    double mag_p2 = sqrt(*sx_p2 * *sx_p2 + *sy_p2 * *sy_p2 + *sz_p2 * *sz_p2);

    struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
    double a1 = o1.a;
    double Om1 = o1.Omega;
    double i1 = o1.inc;
    double pom1 = o1.pomega;
    double e1 = o1.e;

    struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *p2, *sun);
    double a2 = o2.a;
    double Om2 = o2.Omega;
    double i2 = o2.inc;
    double pom2 = o2.pomega;
    double e2 = o2.e;

    fprintf(of, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, sx_sun, sy_sun, sz_sun, mag_sun, sx_p1, sy_p1, sz_p1, mag_p1, sx_p2, sy_p2, sz_p2, mag_p2, a1, Om1, i1, pom1, e1, a2, Om2, i2, pom2, e2);  // prints the spins and orbits of all bodies
    fclose(of);
  }

  if(reb_output_check(sim, 100.*M_PI)){        // outputs to the screen
      reb_output_timing(sim, tmax);
  }
}