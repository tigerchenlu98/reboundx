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

char title[100] = "test2";
char title_stats[100] = "zlk_1112/final_stats_real_Q";
char title_remove[100] = "rm -v test2";
int ind;
int printed_stats=1;

double tmax = 1e5*2*M_PI;
int main(int argc, char* argv[]){
 struct reb_simulation* sim = reb_create_simulation();

 ind = 0;
 if (argc == 2){
   strcat(title, argv[1]);
   strcat(title_remove, argv[1]);
   ind = atoi(argv[1]);
 }

 //sim->rand_seed = ind;

 struct reb_particle star = {0};
 star.m = 0.809;//reb_random_uniform(sim, 0.809 - 0.03, 0.809 + 0.02);
 star.r = 0.683*0.00465;
 /*
 star.x=-1.418804e-02;
 star.y=2.190454e-02;
 star.z=1.837716e-03;
 star.vx=-5.112312e-04;
 star.vy=-8.338613e-04;
 star.vz=5.512050e-04;
 */
 star.x=9.686093e-03;
 star.y=9.954939e-03;
 star.z=-8.402706e-03;
 star.vx=-4.138184e-04;
 star.vy=1.963224e-03;
 star.vz=-3.188257e-04;

 struct reb_particle p1 = {0};
 double planet_m  = 0.0736*9.55e-4;//reb_random_uniform(sim, 0.0736 - 0.0047, 0.0736 + 0.0047) * 9.55e-4;
 double planet_r = 4.36*4.2588e-5;//reb_random_uniform(sim, 4.36 - 0.06, 4.36 + 0.06) * 4.2588e-5;
 p1.m=planet_m;
 p1.r=planet_r;
 /*
 p1.x=3.614459e-01;
 p1.y=-9.446162e-01;
 p1.z=-1.169141e+00;
 p1.vx=-8.264495e-02;
 p1.vy=2.521430e-02;
 p1.vz=-1.269234e-01;
 */
 p1.x=3.203482e-01;
 p1.y=-8.077803e-01;
 p1.z=-9.959403e-01;
 p1.vx=1.816062e-02;
 p1.vy=-9.369794e-02;
 p1.vz=-2.193160e-01;

 struct reb_particle p2 = {0};
 p2.m=3.60 * 9.55e-4;
 /*
 p2.x=3.331215e+00;
 p2.y=-5.135072e+00;
 p2.z=-4.085328e-01;
 p2.vx=1.219881e-01;
 p2.vy=1.957014e-01;
 p2.vz=-1.271099e-01;
 */
 p2.x=-2.285796e+00;
 p2.y=-2.325994e+00;
 p2.z=1.997612e+00;
 p2.vx=9.700483e-02;
 p2.vy=-4.600531e-01;
 p2.vz=7.950706e-02;

 reb_add(sim,star);
 reb_add(sim,p1);
 //reb_add(sim,p2);

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

  //rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default


  struct reb_vec3d star_omega = {0};
  /*
  star_omega.x=-3.656524e-06;
  star_omega.y=-2.859455e-06;
  star_omega.z=1.460000e+01;
  */
  star_omega.x=-8.641459e-07;
  star_omega.y=-6.572553e-06;
  star_omega.z=1.460000e+01;
  rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", star_omega);

  //const double solar_Q = 1e6;
  //double solar_tau = 1 / (2 * solar_Q * orb.n);
  // rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

  // Planet
  rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.1);
  rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * planet_m * planet_r * planet_r);

  struct reb_vec3d planet_omegav = {0};
  /*
  planet_omegav.x=6.494922e+01;
  planet_omegav.y=-3.052679e+02;
  planet_omegav.z=-2.732072e+01;
  */
  planet_omegav.x=1.760493e+03;
  planet_omegav.y=1.032507e+03;
  planet_omegav.z=-3.212885e+02;

  rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", planet_omegav);

  //const double planet_Q = 1000.; //Artificially low
  //struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
  //double planet_n = 2 * np.pi / (0.013388 / (2 * np.pi));
  rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1e-7);


  // add GR precession:
  struct rebx_force* gr = rebx_load_force(rebx, "gr");
  rebx_add_force(rebx, gr);
  rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units

  reb_move_to_com(sim);
  //rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
  rebx_spin_initialize_ode(rebx, effect);

  //system("rm -v test.txt");        // delete previous output file
  system(title_remove);

  FILE* of = fopen(title, "w");
  //fprintf(of, "#Seed: %d,%e,%e,%e,%e,%e,%e\n", index, planet_m, planet_r, planet_a, planet_e, planet_omega, planet_inc);
  //;fprintf(of, "t,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi,mag_p,theta_p,phi_p,n1x,n1y,n1z,n2x,n2y,n2z,sx,sy,sz,Om1,Om2\n");
  fprintf(of, "t,a1,i1,e1,p_ob\n");
  fclose(of);
  tmax=1e6*2*M_PI;
  reb_integrate(sim, tmax);
  rebx_free(rebx);
  reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
   // Output spin and orbital information to file
   if(reb_output_check(sim, 1 * 2 * M_PI)){        // outputs every 100 years
     struct rebx_extras* const rebx = sim->extras;
     FILE* of = fopen(title, "a");
     if (of==NULL){
         reb_error(sim, "Can not open file.");
         return;
     }

     struct reb_particle* sun = &sim->particles[0];
     struct reb_particle* p1 = &sim->particles[1];
     //struct reb_particle* pert = &sim->particles[2];

     // orbits
     struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
     double a1 = o1.a;
     double e1 = o1.e;
     double i1 = o1.inc;
     double Om1 = o1.Omega;
     double pom1 = o1.pomega;
     struct reb_vec3d n1 = o1.hvec;
/*
     struct reb_particle com = reb_get_com_of_pair(sim->particles[0],sim->particles[1]);
     struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *pert, com);
     double a2 = o2.a;
     double e2 = o2.e;
     double i2 = o2.inc;
     double Om2 = o2.Omega;
     double pom2 = o2.pomega;
     struct reb_vec3d n2 = o2.hvec;
*/
     struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

     // Interpret planet spin in the rotating planet frame
     struct reb_vec3d* Omega_p_inv = rebx_get_param(rebx, p1->ap, "Omega");

     // mutual inclination
     double p_ob = obl(*Omega_sun, n1);
     //double pert_ob = obl(*Omega_sun, n2);
     //double mi = obl(n1 , n2);
/*
     if (a1 < 0.05 && printed_stats){
       printed_stats = 0;
       FILE* sf = fopen(title_stats, "a");
       fprintf(sf, "%d,%f,%f,%f\n", ind, e1, p_ob * 180./M_PI, Om1*180./M_PI);
       fclose(sf);
     }
*/
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

     fprintf(of, "%f,%e,%f,%f,%f\n", sim->t,a1,i1,e1,p_ob); // print spins and orbits
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
   if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
       reb_output_timing(sim, tmax);
   }
}
