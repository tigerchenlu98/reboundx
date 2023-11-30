#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "tides_spin.c"

// Kozai for HAT-P-11
void heartbeat(struct reb_simulation* r);

double obl(struct reb_vec3d v1, struct reb_vec3d v2){
  return acos(reb_vec3d_dot(v1,v2) / (sqrt(reb_vec3d_length_squared(v1)) * sqrt(reb_vec3d_length_squared(v2))));
}

char title[100] = "test_";
char title_stats[100] = "zlk_1112/final_stats_real_Q";
char title_remove[100] = "rm -v test_";
int ind;
int printed_stats=1;

double tmax = 1e6*2*M_PI;
int main(int argc, char* argv[]){
 struct reb_simulation* sim = reb_create_simulation();

 ind = 0;
 if (argc == 2){
   strcat(title, argv[1]);
   strcat(title_remove, argv[1]);
   ind = atoi(argv[1]);
 }

 sim->rand_seed = ind;

 struct reb_particle star = {0};
 star.m = 0.809;//reb_random_uniform(sim, 0.809 - 0.03, 0.809 + 0.02);
 star.r = 0.683*0.00465;//reb_random_uniform(sim, 0.683 - 0.009, 0.683 + 0.009) * 0.00465;
 reb_add(sim, star);

 // HAT-P-11b
 // Yee et al 2018
 struct reb_particle planet = {0};
 double planet_m  = reb_random_uniform(sim, 0.0736 - 0.0047, 0.0736 + 0.0047) * 9.55e-4;
 double planet_r = reb_random_uniform(sim, 4.36 - 0.06, 4.36 + 0.06) * 4.2588e-5;
 double planet_a = reb_random_uniform(sim, 0.5, 1.5);
 double planet_e = reb_random_uniform(sim, 0.01, 0.1);
 double planet_Omega = (117.1 - 180.) * (M_PI / 180.); //reb_random_uniform(sim, 0., 2 * M_PI);
 double planet_omega = reb_random_uniform(sim, 0.0, 2 * M_PI);
 double planet_f = reb_random_uniform(sim, 0.0, 2 * M_PI);
 double planet_inc = reb_random_uniform(sim, 6., 39.);

 // HAT-P-11c - treated as a point particle

 struct reb_particle perturber = {0};
 double perturber_a = 4.192;
 double perturber_e = 0.56;
 double perturber_m = 3.06 * 9.55e-4;
 double perturber_inc = 33.5 * M_PI / 180.;
 double perturber_Omega = 117.1 * M_PI / 180.;

 reb_add_fmt(sim, "m r a e inc Omega omega f", planet_m, planet_r, planet_a, planet_e, planet_inc, planet_Omega, planet_omega, planet_f);
 reb_add_fmt(sim, "m a e inc Omega", perturber_m, perturber_a, perturber_e, perturber_inc, perturber_Omega);

 // Reset until over Kozai inclination
 struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
 struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, sim->particles[2], sim->particles[0]);
 double angle = acos(reb_vec3d_dot(o1.hvec, o2.hvec) / (sqrt(reb_vec3d_length_squared(o1.hvec)) * sqrt(reb_vec3d_length_squared(o2.hvec))));
/*
 if (angle < 39. * (180./M_PI)){
   // initial angle is less than the Kozai angle
   // Kill
   printf("Kozai inclination not reached %f %f %f\n", planet_Omega * (180./M_PI), planet_inc * (180./M_PI), angle * (180./M_PI));
   exit(1);
 }
*/
 // Initial conditions
  // Setup constants
  sim->integrator         = REB_INTEGRATOR_IAS15; // IAS15 is used for its adaptive timestep:
  //sim->ri_ias15.adaptive_mode=2;
                                                 // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                 // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
  sim->heartbeat          = heartbeat;

  // Add REBOUNDx effects
  // First tides_spin
  struct rebx_extras* rebx = rebx_attach(sim);

  struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
  rebx_add_force(rebx, effect);


  // Star is roughly Sun-like, shouldn't matter too much
  const double solar_k2 = 0.03;
  rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
  rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.28 * star.m * star.r * star.r);

  const double solar_spin_period = 25. * 2. * M_PI / 365.;
  const double solar_spin = (2 * M_PI) / solar_spin_period;
  rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default

/*
  struct reb_vec3d star_omega = {0};
  star_omega.x=-4.696752e-05;
  star_omega.y=7.965245e-05;
  star_omega.z=1.460000e+01;
  rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", star_omega);
  */
  //const double solar_Q = 1e6;
  //double solar_tau = 1 / (2 * solar_Q * orb.n);
  // rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

  // Planet
  double planet_k2 = reb_random_uniform(sim, 0.1, 0.4);
  rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
  rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * planet_m * planet_r * planet_r);

  const double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
  const double spin_p = (2. * M_PI) / spin_period_p;
  const double theta_p = 0. * M_PI / 180.;
  const double phi_p = 0. * M_PI / 180;
  struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
  rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);
/*
  struct reb_vec3d planet_omegav = {0};
  planet_omegav.x=8.859690e+01;
  planet_omegav.y=4.499880e+02;
  planet_omegav.z=-4.263520e+02;

  rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", planet_omegav);
*/
  //const double planet_Q = 1000.; //Artificially low
  //struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
  //double planet_n = 2 * np.pi / (0.013388 / (2 * np.pi));
  rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1e-4);


  // add GR precession:
  struct rebx_force* gr = rebx_load_force(rebx, "gr");
  rebx_add_force(rebx, gr);
  rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units

  reb_move_to_com(sim);

  // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
  // nodes with the invariable plane (along z cross newz)
  struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
  struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
  struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
  if (isnan(rot.r)) {
    rot = reb_rotation_identity();
  }
  //rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
  rebx_spin_initialize_ode(rebx, effect);

  //system("rm -v test.txt");        // delete previous output file
  system(title_remove);

  FILE* of = fopen(title, "w");
  fprintf(of, "#Seed: %d,%e,%e,%e,%e,%e,%e,%e,%e\n", index, planet_m, planet_r, planet_a, planet_e, planet_omega, planet_inc, planet_f, planet_k2);
  fprintf(of, "t,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi\n");
  //fprintf(of, "t,a1,i1,e1,p_ob,mag_p,theta_p,phi_p\n");
  //"t,ssx,ssy,ssz,mag1,theta1,phi1,a1,e1,nx1,ny1,nz1,nOm1,pom1,a2,e2,i2,Om2,pom2,nx2,ny2,nz2,p_ob,pert_ob,mi\n");
  fclose(of);

  //reb_integrate(sim, 0.41*1e6*2*M_PI);
  //reb_simulationarchive_snapshot(sim, "archive_checked.bin");
  //tmax=509566.;
  reb_integrate(sim, tmax);
  /*
  for (unsigned int i = 0; i < sim->N; i++){
    struct reb_particle* p = &sim->particles[i];
    printf("%d %e %e %e %e %e %e\n", i, p->x,p->y,p->z,p->vx,p->vy,p->vz);
    if (i == 0 || i == 1){
      struct reb_vec3d* Omega = rebx_get_param(rebx, p->ap, "Omega");
      printf("%e %e %e\n", Omega->x, Omega->y, Omega->z);
    }
  }
  */
  //reb_integrate(sim, 1e6*2*M_PI);
  rebx_free(rebx);
  reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
   // Output spin and orbital information to file
   if(reb_output_check(sim, 10. * 2 * M_PI)){        // outputs every 100 years
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

     if (a1 < 0.05258 && printed_stats){
       printed_stats = 0;
       FILE* sf = fopen(title_stats, "a");
       fprintf(sf, "%d,%f,%f\n", ind, e1, p_ob * 180./M_PI);
       fclose(sf);
       exit(1);
     }

     // Transform spin vector into planet frame, w/ z-axis aligned with orbit normal and x-axis aligned with line of nodes
     struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, n1);
     struct reb_rotation rot = reb_rotation_init_to_new_axes(n1, line_of_nodes); // Arguments to this function are the new z and x axes
     if (isnan(rot.r)) {
       rot = reb_rotation_identity();
     }
     struct reb_vec3d srot = reb_vec3d_rotate(*Omega_p_inv, rot); // spin vector in the planet's frame

     // Interpret the spin axis in the more natural spherical coordinates
     /*
     double mag_p;
     double theta_p;
     double phi_p;
     reb_tools_xyz_to_spherical(srot, &mag_p, &theta_p, &phi_p);
     */

     fprintf(of, "%f,%e,%f,%f,%f,%e,%f,%f,%f,%f\n", sim->t,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi); // print spins and orbits
     //fprintf(of, "%f,%e,%f,%f,%f,%e,%f,%f\n", sim->t,a1,i1,e1,p_ob,mag_p,theta_p,phi_p);
     fclose(of);
   }
/*
   if(reb_output_check(sim, 0.01*1e6 * 2 * M_PI)){
     FILE* fe = fopen(title_exact, "a");
     struct rebx_extras* const rebx = sim->extras;
     fprintf(fe, "%f", sim->t);
     for (unsigned int i = 0; i < sim->N; i++){
       struct reb_particle* p = &sim->particles[i];
       fprintf(fe, ",%e,%e,%e,%e,%e,%e", p->x,p->y,p->z,p->vx,p->vy,p->vz);
       if(i==0 || i==1){
         struct reb_vec3d* Omega = rebx_get_param(rebx, p->ap, "Omega");
         fprintf(fe, ",%e,%e,%e",Omega->x,Omega->y,Omega->z);
       }
     }
     fprintf(fe,"\n");
     fclose(fe);

   }
*/
   //if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
  //     reb_output_timing(sim, tmax);
   //}
}
