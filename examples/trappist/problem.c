/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 *
 * This is the same Kozai example implemented in Lu et. al (2023)
 * Also, see the ipython examples prefixed TidesSpin for in-depth exploration of the parameters that can be set in this simulation.
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax = 200 * M_PI * 2.;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Initial conditions
    sim->integrator         = REB_INTEGRATOR_WHFAST; // IAS15 is used for its adaptive timestep:
                                                   // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                   // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
    sim->heartbeat          = heartbeat;

    // Initial conditions
    // Sun
    struct reb_particle star = {0};
    star.m  = 0898;
    star.r = 0.1192 * 0.00465;
    reb_add(sim, star);

    // Planet
    double pm = 10 * 3e-6;
    double pa = 0.5;
    double pr = 3 * 4.26e-5;
    reb_add_fmt(sim, "m a r", pm, pa, pr);

    struct reb_particle planet = sim->particles[1];

    // Test particle
    reb_add_fmt(sim, "primary a", planet, 5. * pr);

    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    sim->dt = orb.P / 10.6789;     // initial timestep as a function of orbital period

    // Add REBOUNDx effects
    // First tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    // Sun
    const double solar_k2 = 0.07;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * star.m * star.r * star.r);

    const double solar_spin_period = 27. * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default

    // Planet
    double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = 2 * M_PI / spin_period_p;
    const double theta_p = 30. * M_PI / 180.;
    const double phi_p = 0. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);

    const double ur_k2 = 0.3;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", ur_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.225 * pm * pr * pr);

    reb_move_to_com(sim);
    rebx_spin_initialize_ode(rebx, effect);

    system("rm -v output.txt");        // delete previous output file
    FILE* of = fopen("output.txt", "w");
    fprintf(of, "t,px,py,pz");
    //for (int i = 1; i < sim->N; i++){
    //  fprintf(of, ",a%d,e%d,i%d,Om%d,pom%d",i,i,i,i,i);
  //  }
    //fprintf(of,"\n");
    fclose(of);

    reb_integrate(sim, tmax);

    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 10)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("output.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }
      struct reb_particle* test = &sim->particles[2];

      fprintf(of, "%e,%e,%e,%e\n", sim->t, test->x, test->y, test->z); // print spins and orbits
      fclose(of);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
