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

//char title[100] = "test_bigr";
char title_stats[100] = "zlk35_Linder_stats";
//char title_remove[100] = "rm -v test_bigr";
int ind;
int printed_stats=1;
double planet_a;
double rfac;
double angle;
double planet_k2;

double tscale = (4./365.) * 2 * M_PI;

double tmax = 1e7*2*M_PI;

// RADIUS INFLATION
double STEFF = 4780.;
double SB = 3.60573e-21;
double SLUMINOSITY;
double FLUX_EARTH;

double get_flux(double r, double l){
    return l / (4. * M_PI * r * r);
}
double R_EARTH = 4.259e-5;
double SUN_TEFF = 5780;

double C0 = 0.131;

double T1_Cs[4] = {-0.348,0.631,0.104,-0.179};

double T2_Cs[4][4] = {{0.209,0.028,-0.168,0.008},
                 {0.,0.086,-0.045,-0.036},
                 {0.,0.,0.052,0.031},
                 {0.,0.,0.,-0.009}};

double get_Rp(double flux){
    double x1 = log10(23.4);
    double x2 = log10(0.2/0.05);
    double x3 = log10(flux/FLUX_EARTH);
    double x4 = log10(5/5);
    double x_array[4] = {x1,x2,x3,x4};

    double t1 = 0.;
    for (int i = 0; i < 4; i++){
        t1 += x_array[i] * T1_Cs[i];
    }

    double t2 = 0;
    for (int i = 0; i < 4; i++){
      for (int j = i; j < 4; j++){
        t2 += T2_Cs[i][j] * x_array[i] * x_array[j];
      }
    }

    double exponent = C0 + t1 + t2;

    return pow(10., exponent) * R_EARTH;
}

// Interpolated curves
// 0.0881, 10 Earth Mass Core
double C1[3] = {2.15926161e-05, 4.02386269e-04, 2.18731715e-03};

// 0.0881, 25 earth mass core
double C2[3] = {3.89300825e-06, 7.32405935e-05, 5.10326839e-04};

// 0.115, 25 earth mass core
double C3[3] = {7.02292131e-06, 1.31862046e-04, 8.59659186e-04};

// Sup-Neptune, icy core, 20 earth masses, 20% envelope fraction
double C4[3] = {1.65036250e-06, 3.70396793e-05, 4.10083278e-04};


double LJUP = 9.43093e-16;
// Linder, interpolation is from Jupiter Luminosities->Earth Radii
double LI[3] = {0.35031788, 1.44509504, 5.63739982};

double interpolate_Rp(double lum){
  if (lum/LJUP < 6e-1){
    return 5.0 * R_EARTH;
  }
  else{
    return (LI[0] * log10(lum/LJUP) * log10(lum/LJUP) + LI[1] * log10(lum/LJUP) + LI[2]) * R_EARTH;
  }
}

double rebx_h1(double e){
  return 1. + (3./2.)*(e*e) + (1./8.)*(e*e*e*e) / pow((1. - e*e), (9./2.));
}

double rebx_h2(double e){
  return 1. + (9./2.)*(e*e) + (5./8.)*(e*e*e*e) / pow((1. - e*e), (9./2.));
}

double rebx_h3(double e){
  return ((1. + 3. * (e*e) + 3. * (e*e*e*e))/8.) / pow((1. - e*e), (9./2.));
}

double rebx_h4(double e){
  return (1. + (15. * (e*e)/2.) + (45. * (e*e*e*e)/8.) * (5. * (e*e*e*e*e*e)/16.)) / pow((1. - e*e), 6.);
}

double rebx_h5(double e){
  return (1. + (31. * (e*e)/2.) + (255. * (e*e*e*e)/8.) + (185. * (e*e*e*e*e*e)/16.) + (25. * (e*e*e*e*e*e*e*e)/64.)) / pow((1. - e*e), 15./2.);
}

double Etide(struct reb_simulation* sim, struct reb_extras* rebx, struct reb_particle* planet, struct reb_particle* star){
  struct reb_orbit o = reb_orbit_from_particle(sim->G, *planet, *star);
  const double m1 = star->m;
  const double m2 = planet->m;
  const double r2 = planet->r;
  const double a = o.a;
  const double e = o.e;

  struct reb_vec3d hhat = reb_vec3d_normalize(o.hvec);
  struct reb_vec3d ehat = reb_vec3d_normalize(o.evec);
  struct reb_vec3d qhat = reb_vec3d_cross(hhat, ehat);

  const double* k2 = rebx_get_param(rebx, planet->ap, "k2");
  //const double* tau = rebx_get_param(rebx, planet->ap, "tau");
  const double tau = 1e-8; // manually set this
  const double invQ = (2. * o.n * tau);
  struct reb_vec3d* Omega = rebx_get_param(rebx, planet->ap, "Omega");

  const double prefactor = (m1*m2)/(m1+m2) * a*a * o.n * (m1/m2) * pow((r2/a), 5.) * 6. * (*k2) * invQ;

  const double Omega_e = reb_vec3d_dot(*Omega, ehat);
  const double Omega_q = reb_vec3d_dot(*Omega, qhat);
  const double Omega_h = reb_vec3d_dot(*Omega, hhat);

  const double t1 = 0.5 * (Omega_e * Omega_e * rebx_h1(e) + Omega_q * Omega_q * rebx_h2(e));
  const double t2 = Omega_h * Omega_h * rebx_h3(e);
  const double t3 = -2. * o.n * Omega_h * rebx_h4(e);
  const double t4 = o.n*o.n + rebx_h5(e);

  return fabs(prefactor * (t1+t2+t3+t4));
}

int main(int argc, char* argv[]){
 struct reb_simulation* sim = reb_simulation_create();

 //FLUX_EARTH = SB * 4 * M_PI * 0.00465 * 0.004652 * pow(SUN_TEFF, 4.) / (4. * M_PI);

 ind = 0;
 if (argc == 2){
   //strcat(title, argv[1]);
   //strcat(title_remove, argv[1]);
   ind = atoi(argv[1]);
 }

 sim->rand_seed = ind;

 struct reb_particle star = {0};
 star.m = 0.809;//reb_random_uniform(sim, 0.809 - 0.03, 0.809 + 0.02);
 star.r = 0.683*0.00465;//reb_random_uniform(sim, 0.683 - 0.009, 0.683 + 0.009) * 0.00465;

 //double star_area = 4 * M_PI * star.r * star.r;
 //SLUMINOSITY = SB * star_area * pow(STEFF, 4.);
 reb_simulation_add(sim, star);

 // HAT-P-11b
 // Yee et al 2018
 struct reb_particle planet = {0};
 double planet_m  = reb_random_uniform(sim, 0.0736 - 0.0047, 0.0736 + 0.0047) * 9.55e-4;
 //rfac = reb_random_uniform(sim, 1.2, 3.0);
 double planet_r = 4.36 * 4.2588e-5;// doesn't matter, instantly overwritten
 planet_a = reb_random_uniform(sim, 0.3, 0.8);
 double planet_e = reb_random_uniform(sim, 0.01, 0.1);
 double planet_Omega = (117.1 - 180.) * (M_PI / 180.); //reb_random_uniform(sim, 0., 2 * M_PI);
 double planet_omega = 0.;//reb_random_uniform(sim, 0.0, 2 * M_PI);
 double planet_f = 0;//reb_random_uniform(sim, 0.0, 2 * M_PI);
 double planet_inc = reb_random_uniform(sim, 6. * M_PI/180., 39. * M_PI/180.);
 //double planet_inc = 106. * M_PI/180.;

 // HAT-P-11c - treated as a point particle

 //struct reb_particle perturber = {0};
 double perturber_a = 4.192;//reb_random_uniform(sim, 4.192 - 0.072, 4.192 + 0.072);
 double perturber_e = 0.56;//reb_random_uniform(sim, 0.56 - 0.035, 0.56 + 0.035);
 double perturber_m = 3.06 * 9.55e-4;//reb_random_uniform(sim, 3.06 - 0.41, 3.06 + 0.43) * 9.55e-4;
 double perturber_inc = 33.5 * M_PI/180.;//reb_random_uniform(sim, 33.5 - 4.4, 33.5 + 6.1) * M_PI / 180.;
 double perturber_Omega = 117.1 * M_PI / 180.;

 reb_simulation_add_fmt(sim, "m r a e inc Omega omega f", planet_m, planet_r, planet_a, planet_e, planet_inc, planet_Omega, planet_omega, planet_f);
 reb_simulation_add_fmt(sim, "m a e inc Omega", perturber_m, perturber_a, perturber_e, perturber_inc, perturber_Omega);

 // Reset until over Kozai inclination
   struct reb_orbit o1 = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
   struct reb_orbit o2 = reb_orbit_from_particle(sim->G, sim->particles[2], sim->particles[0]);
   angle = acos(reb_vec3d_dot(o1.hvec, o2.hvec) / (sqrt(reb_vec3d_length_squared(o1.hvec)) * sqrt(reb_vec3d_length_squared(o2.hvec))));
   //printf("%f\n", angle * 180./M_PI);

 // Initial conditions
  // Setup constants
  sim->integrator         = REB_INTEGRATOR_IAS15; // IAS15 is used for its adaptive timestep:
  sim->ri_ias15.adaptive_mode = 2;
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
  rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", 1e-8);

  // Planet
  planet_k2 = reb_random_uniform(sim, 0.1, 0.8);
  rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
  rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * planet_m * planet_r * planet_r);

  const double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
  const double spin_p = (2. * M_PI) / spin_period_p;
  struct reb_orbit ob = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
  struct reb_vec3d Omega_sv = reb_vec3d_mul(reb_vec3d_normalize(ob.hvec), spin_p);
  rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);
  rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1e-5);


  // add GR precession:
  struct rebx_force* gr = rebx_load_force(rebx, "gr_potential");
  rebx_add_force(rebx, gr);
  rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units
  reb_simulation_move_to_com(sim);

  // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
  // nodes with the invariable plane (along z cross newz)
  struct reb_vec3d newz = reb_vec3d_add(reb_simulation_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
  struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
  struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
  if (isnan(rot.r)) {
    rot = reb_rotation_identity();
  }
  //rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
  rebx_spin_initialize_ode(rebx, effect);

  //system("rm -v test.txt");        // delete previous output file
  double lum = Etide(sim, rebx, &sim->particles[1], &sim->particles[0]);
  double re = interpolate_Rp(lum);
  //printf("Etid: %f %f\n", lum/LJUP, re/R_EARTH);
  //exit(1);
  //system(title_remove);
/*
  FILE* of = fopen(title, "w");
  fprintf(of, "#Seed: %d,%e,%e,%e,%e,%e,%e,%e,%e\n", index, planet_m, planet_r, planet_a, planet_e, planet_omega, planet_inc, planet_f, planet_k2);
  fprintf(of, "t,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi,mag_p,theta_p,pr\n");
  //fprintf(of, "t,a1,i1,e1,p_ob,mag_p,theta_p,phi_p\n");
  //"t,ssx,ssy,ssz,mag1,theta1,phi1,a1,e1,nx1,ny1,nz1,nOm1,pom1,a2,e2,i2,Om2,pom2,nx2,ny2,nz2,p_ob,pert_ob,mi\n");
  fclose(of);
*/
  reb_simulation_integrate(sim, tmax);
  struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
  struct reb_vec3d n1 = orb.hvec;

  struct reb_particle* sun = &sim->particles[0];
  struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

  double p_ob = obl(*Omega_sun, n1);

  FILE* sf = fopen(title_stats, "a");
  fprintf(sf, "%d,%f,%f,%f,%f,%f,%f,1\n", ind, orb.e, sim->particles[1].r / R_EARTH, planet_a, planet_k2, angle*180./M_PI,p_ob * 180./M_PI);
  fclose(sf);
  rebx_free(rebx);
  reb_simulation_free(sim);

}

void heartbeat(struct reb_simulation* sim){
   // Radius INFLATION
   if(reb_simulation_output_check(sim, tscale)){
     struct rebx_extras* const rebx = sim->extras;
     struct reb_particle* sun = &sim->particles[0];
     struct reb_particle* p1 = &sim->particles[1];
  /*
     double dx = p1->x - sun->x;
     double dy = p1->y - sun->y;
     double dz = p1->z - sun->z;
     double d = sqrt(dx*dx+dy*dy+dz*dz);

     double flux = get_flux(d, SLUMINOSITY);
     */
     double luminosity = Etide(sim, rebx, p1, sun);
     double rad = interpolate_Rp(luminosity);
     p1->r = rad;


     // Check for break
     struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
     struct reb_vec3d n1 = orb.hvec;

     //struct reb_particle* sun = &sim->particles[0];
     struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

     double p_ob = obl(*Omega_sun, n1);
     if (orb.a < 0.05){
       FILE* sf = fopen(title_stats, "a");
       fprintf(sf, "%d,%f,%f,%f,%f,%f,%f,0\n", ind, orb.e, sim->particles[1].r / R_EARTH, planet_a, planet_k2, angle*180./M_PI,p_ob * 180./M_PI);
       fclose(sf);
       exit(1);
     }

   }
   // Output spin and orbital information to file
/*
   if(reb_simulation_output_check(sim, 10. * 2 * M_PI)){        // outputs every 100 years
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

     //struct reb_particle com = reb_get_com_of_pair(sim->particles[0],sim->particles[1]);
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
     double mi = obl(n1 , n2);

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


     FILE* of = fopen(title, "a");
     fprintf(of, "%f,%e,%f,%f,%f,%e,%f,%f,%f,%f,%e\n", sim->t,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi,mag_p,theta_p,p1->r); // print spins and orbits
     //fprintf(of, "%f,%e,%f,%f,%f,%e,%f,%f\n", sim->t,a1,i1,e1,p_ob,mag_p,theta_p,phi_p);
     fclose(of);

   }


   if(reb_simulation_output_check(sim, 20.*M_PI)){        // outputs to the screen
       reb_simulation_output_timing(sim, tmax);
   }
   */

}
