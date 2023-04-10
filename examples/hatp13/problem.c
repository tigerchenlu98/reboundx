 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <math.h>
 #include <string.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax;
double k2;
int i_glob;
int j_glob;
char title[100] = "out.txt";

double obl(struct reb_vec3d v1, struct reb_vec3d v2){
  return acos(reb_vec3d_dot(v1,v2) / (sqrt(reb_vec3d_length_squared(v1)) * sqrt(reb_vec3d_length_squared(v2))));
}

int main(int argc, char* argv[]){
  double k2s[10] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
  double mincs[21] = {-50., -45., -40., -35., -30., -25., -20., -15., -10., -5., 0.0, 5., 10., 15., 20., 25., 30., 35., 40., 45., 50.};

  system("rm -v out.txt");        // delete previous output file
  FILE* of = fopen(title, "w");

  struct reb_particle star = {0};
  star.m  = 1.22;
  star.r = 1.56 * 0.00465;

  double planet_m  = 0.851 * 9.55e-4; // A Jupiter-like planet
  double planet_r = 1.28 * 4.676e-4;
  double planet_a = 0.0426;
  double planet_e = 0.04;//0.021;
  double planet_inc = 0.0 * M_PI / 180.;//0.021;

  // The perturber - treated as a point particle
  double perturber_a = 1.186;
  double perturber_e = 0.691;
  double perturber_m;
  double minc;
  double perturber_inc;

  //double k2s[3] =  {0.1, 0.2, 0.3};
  //double mincs[2] = {0.0, 10.};

  for (int i = 0; i < 10; i++){
    i_glob = i;
    for (int j = 0; j < 21; j++){
        j_glob = j;
        struct reb_simulation* sim = reb_create_simulation();
        // Initial conditions
        // Setup constants
        sim->integrator         = REB_INTEGRATOR_WHFAST; // IAS15 is used for its adaptive timestep:
                                                       // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                       // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
        sim->heartbeat          = heartbeat;

        // Initial conditions
        reb_add(sim, star);

        // struct reb_particle planet = {0};
        reb_add_fmt(sim, "m r a e inc", planet_m, planet_r, planet_a, planet_e, planet_inc);

        minc = mincs[j];
        printf("%f\n", minc);
        perturber_inc = minc * M_PI / 180.;
        perturber_m  = 15.2 * 9.55e-4 / cos(fabs(planet_inc - perturber_inc));
        reb_add_fmt(sim, "m a e inc", perturber_m, perturber_a, perturber_e, perturber_inc);

        struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
        sim->dt = orb.P / 10.12345;     // initial timestep as a function of orbital period

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
        // rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

        // Planet
        const double planet_k2 = k2s[i];
        k2 = planet_k2;
        rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
        rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.51 * planet_m * planet_r * planet_r);

        const double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
        const double spin_p = (2. * M_PI) / spin_period_p;
        const double theta_p = 0. * M_PI / 180.;
        const double phi_p = 0. * M_PI / 180;
        struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
        rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);

        const double planet_Q = 10;
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
        if isnan(rot.r) {
          rot = reb_rotation_identity();
        }
        rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)

        rebx_spin_initialize_ode(rebx, effect);
        //fprintf(of, "t,ssx,ssy,ssz,mag1,theta1,phi1,a1,e1,nx1,ny1,nz1,nOm1,pom1,a2,e2,i2,Om2,pom2,nx2,ny2,nz2,p_ob,pert_ob,mi\n");

        // Simulation Archive
        //system("rm -v archive.bin");
        //reb_simulationarchive_automate_step(sim, "archive.bin", 20000);

        tmax = 6. * ((2. / (21 * orb.n)) * (planet_Q / planet_k2) * (planet_m / star.m) * pow((orb.a / planet_r), 5));
        //printf("tmax = %f years\n", tmax / 2 * M_PI);

        reb_integrate(sim, tmax);

        rebx_free(rebx);
        reb_free_simulation(sim);
      }
    }
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 100)){        // outputs every 100 years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen(title, "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* p1 = &sim->particles[1];
      struct reb_particle* pert = &sim->particles[2];

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

      struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

      // Interpret planet spin in the rotating planet frame
      struct reb_vec3d* Omega_p_inv = rebx_get_param(rebx, p1->ap, "Omega");

      // mutual inclination
      double p_ob = obl(*Omega_sun, n1);
      double pert_ob = obl(*Omega_sun, n2);
      double mi = obl(n1 , n2);

      // Transform spin vector into planet frame, w/ z-axis aligned with orbit normal and x-axis aligned with line of nodes
      struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, n1);
      struct reb_rotation rot = reb_rotation_init_to_new_axes(n1, line_of_nodes); // Arguments to this function are the new z and x axes
      if isnan(rot.r) {
        rot = reb_rotation_identity();
      }
      struct reb_vec3d srot = reb_vec3d_rotate(*Omega_p_inv, rot); // spin vector in the planet's frame

      // Interpret the spin axis in the more natural spherical coordinates
      double mag_p;
      double theta_p;
      double phi_p;
      reb_tools_xyz_to_spherical(srot, &mag_p, &theta_p, &phi_p);

      fprintf(of, "%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", i_glob,j_glob,sim->t, Omega_sun->x, Omega_sun->y, Omega_sun->z, mag_p, theta_p, phi_p, a1, e1, i1, n1.x, n1.y, n1.z, Om1, pom1, a2, e2, i2, n2.x, n2.y, n2.z, Om2, pom2, p_ob, pert_ob, mi, k2); // print spins and orbits

      fclose(of);
    }

    //if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
    //    reb_output_timing(sim, tmax);
    //}
}
