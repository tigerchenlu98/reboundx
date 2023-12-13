/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 *
 * This example includes self-consistent spin, tidal & dynamical effects
 * as well as general relativity
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <math.h>
 #include <time.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax = 1e7 * 2 * M_PI;
clock_t begin;
clock_t end;
double e_init;
char title[100] = "ias15_tides";
char title_remove[100] = "rm -rf ias15_tides";

double obl(struct reb_vec3d v1, struct reb_vec3d v2){
  return acos(reb_vec3d_dot(v1,v2) / (sqrt(reb_vec3d_length_squared(v1)) * sqrt(reb_vec3d_length_squared(v2))));
}

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();
    // Setup constants
    //sim->dt             = 2*M_PI*1.;     // initial timestep
    sim->integrator        = REB_INTEGRATOR_TRACE;
    //sim->ri_ias15.adaptive_mode = 2;
    sim->heartbeat        = heartbeat;

    // Initial conditions

    struct reb_particle star = {0};
    star.m  = 1.1;
    star.r = 0.00465;
    reb_simulation_add(sim, star);

    double planet_m  = 7.8 * 9.55e-4; // in Jupiter masses
    double planet_r = 4.676e-4;
    double planet_a = 5.0;
    double planet_e = 0.1;
    reb_simulation_add_fmt(sim, "m r a e", planet_m, planet_r, planet_a, planet_e);

    struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    sim->dt = o.P/20.12345;
    sim->exact_finish_time=0;

    // The perturber
    double perturber_inc = 85.6 * (M_PI / 180.);
    double perturber_mass = 1.1;
    double perturber_a  = 1000.;
    double perturber_e = 0.5;
    reb_simulation_add_fmt(sim, "m a e inc", perturber_mass, perturber_a, perturber_e, perturber_inc);

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    // Sun
    const double solar_spin_period = 20 * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_k2 = 0.028;
    //const double solar_tau = (0.2 / solar_k2) * (1 / 86400.) * 2 * M_PI / 365.; // seconds to reb years
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.08 * star.m * star.r * star.r);
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin});

    double solar_Q = 1e6;
    struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    double solar_tau = 1. / (2. * solar_Q * orb.n);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

    // P1
    const double spin_period_p = (10./24.) * 2. * M_PI / 365.; // days to reb years
    const double spin_p = (2. * M_PI) / spin_period_p;
    const double planet_k2 = 0.51;
    const double planet_Q = 3. * 1e5;
    // const double planet_tau = (0.02 / planet_k2) * (1 / 86400.) * 2 * M_PI / 365.; // seconds to reb years
    const double theta_1 = 0. * M_PI / 180.;
    const double phi_1 = 0. * M_PI / 180;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * planet_m * planet_r * planet_r);

    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_1, phi_1);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);
    double planet_tau = 1. / (2. * planet_Q * orb.n);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", planet_tau);

    // add GR:
    struct rebx_force* gr = rebx_load_force(rebx, "gr");
    rebx_add_force(rebx, gr);

    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units
    reb_simulation_move_to_com(sim);

    struct reb_vec3d newz = reb_vec3d_add(reb_simulation_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    //struct reb_vec3d newz = reb_simulation_angular_momentum(sim);
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
    //reb_simulation_irotate(sim, rot);

    rebx_spin_initialize_ode(rebx, effect);

    system(title_remove);
    FILE* f = fopen(title,"w");
    fprintf(f, "t,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi,theta_p,time_spent\n");
    //fprintf(f, "t,a1,i1,e1,a2,i2,e2,mi,time_spent\n");

    begin = clock();
    e_init = reb_simulation_energy(sim);

    while (sim->t < tmax){
      reb_simulation_integrate(sim, sim->t + 1e5*2*M_PI);
      end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

      struct rebx_extras* const rebx = sim->extras;
      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* p1 = &sim->particles[1];
      struct reb_particle* pert = &sim->particles[2];

      // orbits
      struct reb_orbit o1 = reb_orbit_from_particle(sim->G, *p1, *sun);
      double a1 = o1.a;
      double e1 = o1.e;
      double i1 = o1.inc;
      double Om1 = o1.Omega;
      double pom1 = o1.pomega;
      struct reb_vec3d n1 = o1.hvec;

      struct reb_orbit o2 = reb_orbit_from_particle(sim->G, *pert, *sun);
      double a2 = o2.a;
      double e2 = o2.e;
      double i2 = o2.inc;
      double Om2 = o2.Omega;
      double pom2 = o2.pomega;
      struct reb_vec3d n2 = o2.hvec;
      struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

      // Interpret planet spin in the rotating planet frame
      struct reb_vec3d* Omega_p_inv = rebx_get_param(rebx, p1->ap, "Omega");

      // mutual inclination
      double p_ob = obl(*Omega_sun, n1);
      double pert_ob = obl(*Omega_sun, n2);
      double mi = obl(n1,n2);

      // Transform spin vector into planet frame, w/ z-axis aligned with orbit normal and x-axis aligned with line of nodes
      struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, n1);
      struct reb_rotation rot = reb_rotation_init_to_new_axes(n1, line_of_nodes); // Arguments to this function are the new z and x axes
      if (isnan(rot.r)) {
        rot = reb_rotation_identity();
      }
      struct reb_vec3d srot = reb_vec3d_rotate(*Omega_p_inv, rot); // spin vector in the planet's frame

      // Interpret the spin axis in the more natural spherical coordinates
      double mag_p;
      double theta_p;
      double phi_p;
      reb_tools_xyz_to_spherical(srot, &mag_p, &theta_p, &phi_p);

      fprintf(f, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", sim->t,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi,theta_p,time_spent); // print spins and orbits
      //fprintf(of, "%f,%e,%f,%f,%f,%f,%f,%f,%f,%f\n", sim->t,fabs((reb_simulation_energy(sim) - e_init) / e_init),a1,i1,e1,a2,i2,e2,mi,time_spent); // print spins and orbits
      }
    fclose(f);
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%f\n", time_spent);
    rebx_free(rebx);
    reb_simulation_free(sim);
}

void heartbeat(struct reb_simulation* sim){
   // Output spin and orbital information to file

   if (reb_simulation_output_check(sim, 1e3*2.*M_PI)){
        reb_simulation_output_timing(sim, tmax);
   }
}
