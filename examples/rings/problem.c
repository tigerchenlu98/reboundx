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
double tmax = 1000 * M_PI * 2.;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Initial conditions
    // Setup constants
    sim->integrator         = REB_INTEGRATOR_IAS15; // IAS15 is used for its adaptive timestep:
                                                   // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                   // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
    sim->heartbeat          = heartbeat;
    sim->N_active = 2;

    // Initial conditions
    // Sun
    struct reb_particle star = {0};
    star.m  = 1.0;
    star.r = 1.0 * 0.00465;
    reb_add(sim, star);

    // Planet
    double pm = 10 * 3e-6;
    double pa = 0.5;
    double pr = 3 * 4.26e-5;
    double pinc = 0. * (M_PI / 180.);
    reb_add_fmt(sim, "m a r inc", pm, pa, pr, pinc);

    struct reb_particle planet = sim->particles[1];
    // Test particle
    double tinc = 5. * (M_PI / 180.);
    reb_add_fmt(sim, "primary a inc", planet, 5. * pr, tinc);
    reb_add_fmt(sim, "primary a inc", planet, 30. * pr, tinc);


    // Add REBOUNDx effects
    // First tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    // Sun
    const double solar_k2 = 0.01;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * star.m * star.r * star.r);

    const double solar_spin_period = 27. * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default

    // Planet

    double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = 2 * M_PI / spin_period_p;
    const double theta_p = 20. * M_PI / 180.;
    const double phi_p = 0. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);

    const double p_k2 = 0.7;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", p_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.225 * pm * pr * pr);


    reb_move_to_com(sim);
    rebx_spin_initialize_ode(rebx, effect);


    system("rm -v output.txt");        // delete previous output file
    FILE* of = fopen("output.txt", "w");
    fprintf(of, "t,sx,sy,sz,nx,ny,nz,nt1x,nt1y,nt1z,nt2x,nt2y,nt2z\n");
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
    if(reb_output_check(sim, 100)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("output.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }
      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* p = &sim->particles[1];
      struct reb_particle* t1 = &sim->particles[2];
      struct reb_particle* t2 = &sim->particles[3];

      struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *p, *sun);
      struct reb_vec3d n = o.hvec;

      struct reb_orbit ot1 = reb_tools_particle_to_orbit(sim->G, *t1, *p);
      struct reb_vec3d nt1 = ot1.hvec;

      struct reb_orbit ot2 = reb_tools_particle_to_orbit(sim->G, *t2, *p);
      struct reb_vec3d nt2 = ot2.hvec;

      struct reb_vec3d* s = rebx_get_param(rebx, p->ap, "Omega");
      //struct reb_vec3d* s = {0};

      fprintf(of, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t,s->x,s->y,s->z,n.x,n.y,n.z,nt1.x,nt1.y,nt1.z,nt2.x,nt2.y,nt2.z); // print spins and orbits
      //fprintf(of, "%e,%e,%e,%e,%e,%e,%e,%f\n", sim->t,0.0,0.0,0.0,n.x,n.y,n.z,o.a); // print spins and orbits
      fclose(of);
    }

    //if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
    //    reb_output_timing(sim, tmax);
    //}
}
