 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax = 1e5 * 2 * M_PI; // kept short to run quickly.
                   // set to 3e5 for a full cycle
                   // or 7e6 * 2 * M_PI to reproduce the paper plot

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Initial conditions
    // Setup constants
    sim->integrator         = REB_INTEGRATOR_IAS15;
    sim->heartbeat          = heartbeat;

    // Initial conditions
    // All units are in default REBOUND units. So Solar Masses, AU, and yr/2 pi.
    struct reb_particle star = {0};
    star.m  = 1.21;
    star.r = 1.56 * 0.00465;
    reb_add(sim, star);

    // struct reb_particle planet = {0};
    double planet_m  = 0.85 * 9.55e-4; // A Jupiter-like planet
    double planet_r = 1.5 * 4.676e-4;
    double planet_a = 0.043;
    double planet_e = 0.005;
    double planet_inc = 0.0 * (M_PI / 180.);
    double planet_pomega = 0. * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc pomega", planet_m, planet_r, planet_a, planet_e, planet_inc, planet_pomega);

    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    //sim->dt = orb.P / 15.12345;

    // The perturber - treated as a point particle
    double perturber_m  = 15.2 * 9.55e-4;
    double perturber_a = 1.7;
    double perturber_e = 0.1;
    double perturber_inc = 15.0 * (M_PI / 180.);
    //double perturber_inc = 0.5 * (M_PI / 180.);
    reb_add_fmt(sim, "m a e inc", perturber_m, perturber_a, perturber_e, perturber_inc);

    // The scattered body - treated as a point particle
    double eject_m  = 12. * 9.55e-4;
    double eject_a = 3.097;
    double eject_e = 0.05;
    double eject_inc = 15.0 * (M_PI / 180.);
    //double perturber_inc = 0.5 * (M_PI / 180.);
    reb_add_fmt(sim, "m a e inc", eject_m, eject_a, eject_e, eject_inc);

    // Add REBOUNDx effects
    // First tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);


    // Sun
    const double solar_k2 = 0.03;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.28 * star.m * star.r * star.r);

    const double solar_spin_period = 25. * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default

    const double solar_Q = 1e6;
    // In the case of a spin that is synchronous with a circular orbit, tau is related to the tidal quality factor Q through the orbital mean motion n (see Lu et al. 2023 for discussion). Clearly that's not the case here, but gives us a reasonable starting point to start turning this knob
    double solar_tau = 1 / (2 * solar_Q * orb.n);
    //rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

    // Planet
    const double planet_k2 = 0.25;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.51 * planet_m * planet_r * planet_r);

    const double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = (2. * M_PI) / spin_period_p;

    struct reb_vec3d norm = orb.hvec;
    struct reb_vec3d nhat = reb_vec3d_normalize(norm);
    // struct reb_vec3d Omega_p = reb_vec3d_mul(nhat, spin_rate);
    struct reb_vec3d Omega_p = reb_vec3d_mul(nhat, spin_p); // begin with 0 obliquity
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_p);

    const double planet_Q = 40.;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2.*planet_Q*orb.n));


    // add GR precession:
    struct rebx_force* gr = rebx_load_force(rebx, "gr_potential");
    rebx_add_force(rebx, gr);
    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units

    reb_move_to_com(sim);

    // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
    // nodes with the invariable plane (along z cross newz)
    struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    //rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)

    rebx_spin_initialize_ode(rebx, effect);

    //system("rm -v output.txt");        // delete previous output file
    FILE* of = fopen("scattering.txt", "a");
    fprintf(of, "t,ssx,ssy,ssz,sbx,sby,sbz,a1,e1,i1,Om1,pom1,nbx,nby,nbz,a2,e2,i2,Om2,pom2,ncx,ncy,ncz,ad,ed\n");
    fclose(of);

    reb_integrate(sim, tmax);

    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 10)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("scattering.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* p1 = &sim->particles[1];
      struct reb_particle* pert = &sim->particles[2];
      struct reb_particle* ej = &sim->particles[3];

      // orbits
      struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
      double a1 = o1.a;
      double e1 = o1.e;
      double i1 = o1.inc;
      double Om1 = o1.Omega;
      double pom1 = o1.pomega;
      struct reb_vec3d n1 = o1.hvec;

      struct reb_particle com = reb_get_com_of_pair(sim->particles[0],sim->particles[1]);
      struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *pert, com);
      double a2 = o2.a;
      double e2 = o2.e;
      double i2 = o2.inc;
      double Om2 = o2.Omega;
      double pom2 = o2.pomega;
      struct reb_vec3d n2 = o2.hvec;

      struct reb_orbit o3 = reb_tools_particle_to_orbit(sim->G, *ej, com);
      double a3 = o3.a;
      double e3 = o3.e;

      struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

      // Interpret planet spin in the rotating planet frame
      struct reb_vec3d* Omega_p_inv = rebx_get_param(rebx, p1->ap, "Omega");

      // Transform spin vector into planet frame, w/ z-axis aligned with orbit normal and x-axis aligned with line of nodes
      //struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, n1);
      //struct reb_rotation rot = reb_rotation_init_to_new_axes(n1, line_of_nodes); // Arguments to this function are the new z and x axes
      //struct reb_vec3d srot = reb_vec3d_rotate(*Omega_p_inv, rot); // spin vector in the planet's frame

      // Interpret the spin axis in the more natural spherical coordinates
      //double mag_p;
      //double theta_p;
      //double phi_p;
      //reb_tools_xyz_to_spherical(srot, &mag_p, &theta_p, &phi_p);

      fprintf(of, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, Omega_sun->x, Omega_sun->y, Omega_sun->z, Omega_p_inv->x, Omega_p_inv->y, Omega_p_inv->z, a1, e1, i1, Om1, pom1,n1.x,n1.y,n1.z,a2, e2, i2, Om2, pom2, n2.x,n2.y,n2.z, a3, e3); // print spins and orbits

      fclose(of);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
