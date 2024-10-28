#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "tides_spin.c"

// ZLK for HAT-P-11. Used to generate Figures 7 & 8 in Lu+ (2024)
void heartbeat(struct reb_simulation* r);

double obl(struct reb_vec3d v1, struct reb_vec3d v2){
  return acos(reb_vec3d_dot(v1,v2) / (sqrt(reb_vec3d_length_squared(v1)) * sqrt(reb_vec3d_length_squared(v2))));
}

double planet_k2 = 0.5;
double QPRIME = 3e5;
double tscale = (4./365.) * 2 * M_PI;
double tmax = 1e8*2*M_PI;
double R_EARTH = 4.259e-5;
int ind;

// RADIUS INFLATION
// Calculate tidal luminosity
double rebx_h1(double e){
  return 1. + (3./2.)*(e*e) + (1./8.)*(e*e*e*e) / pow((1. - e*e), (9./2.));
}

double rebx_h2(double e){
  return 1. + (9./2.)*(e*e) + (5./8.)*(e*e*e*e) / pow((1. - e*e), (9./2.));
}

double rebx_h3(double e){
  return ((1. + 3. * (e*e) + 3. * (e*e*e*e)/8.) / pow((1. - e*e), (9./2.)));
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

  struct reb_vec3d* Omega = rebx_get_param(rebx, planet->ap, "Omega");

  const double prefactor = (m1*m2)/(m1+m2) * a*a * o.n * (m1/m2) * pow((r2/a), 5.) * 3./(2. * QPRIME);

  const double Omega_e = reb_vec3d_dot(*Omega, ehat);
  const double Omega_q = reb_vec3d_dot(*Omega, qhat);
  const double Omega_h = reb_vec3d_dot(*Omega, hhat);

  const double t1 = 0.5 * (Omega_e * Omega_e * rebx_h1(e) + Omega_q * Omega_q * rebx_h2(e));
  const double t2 = Omega_h * Omega_h * rebx_h3(e);
  const double t3 = -2. * o.n * Omega_h * rebx_h4(e);
  const double t4 = o.n*o.n * rebx_h5(e);

  return fabs(prefactor * (t1+t2+t3+t4));
}

// MPB interpolation
double MPB[5] = {5.23017622e-04, 3.36801133e-02, 8.00957870e-01, 8.35709522e+00, 3.59588034e+01};
double LSOLAR =  1.08827e-6;

// Fit REBOUND luminosities to earth radii
double interpolate_mpb(double lum){
  double llum = log10(lum/LSOLAR);
  return (MPB[0] * llum * llum * llum * llum + MPB[1] * llum * llum * llum + MPB[2] * llum * llum + MPB[3] * llum + MPB[4]) * R_EARTH;
}

// Pseudo synchronous spin rate
double pseudo_sync(double e, double n){
    return ((1. + 15./2. * e*e + 45./8. * e*e*e*e + 5./16. * e*e*e*e*e*e) / ((1 + 3*e*e + 3/8*e*e*e*e)*pow(1-e*e, 3./2))*n);
}

int main(int argc, char* argv[]){
 struct reb_simulation* sim = reb_simulation_create();

 ind = 0;
 if (argc == 2){
   ind = atoi(argv[1]);
 }

 sim->rand_seed = ind;

 struct reb_particle star = {0};
 star.m = 0.809;
 star.r = 0.683*0.00465;

 reb_simulation_add(sim, star);

 // HAT-P-11b
 // Yee et al 2018
 struct reb_particle planet = {0};
 double planet_m  = 0.0736 * 9.55e-4;
 double planet_r = 4.36 * 4.2588e-5;
 planet_a = reb_random_uniform(sim, 0.154, 0.518);
 double planet_e = reb_random_uniform(sim, 0.01, 0.7);
 double planet_Omega = reb_random_uniform(sim, 0.0, 2 * M_PI);
 double planet_omega = reb_random_uniform(sim, 0.0, 2 * M_PI);
 double planet_f = reb_random_uniform(sim, 0.0, 2 * M_PI);
 double planet_inc = reb_random_uniform(sim, 1.01 * M_PI/180., 97.0 * M_PI/180.);

 // HAT-P-11c - treated as a point particle
 double perturber_a = 4.10;
 double perturber_e = 0.652;
 double perturber_m = 2.68 * 9.55e-4;
 double perturber_inc = reb_random_uniform(sim, 33.3, 50.26) * M_PI / 180.;
 double perturber_Omega = reb_random_uniform(sim, 0.0, 2*M_PI);
 double perturber_omega = reb_random_uniform(sim, 0.0, 2 * M_PI);
 double perturber_f = reb_random_uniform(sim, 0.0, 2 * M_PI);

 reb_simulation_add_fmt(sim, "m r a e inc Omega omega f", planet_m, planet_r, planet_a, planet_e, planet_inc, planet_Omega, planet_omega, planet_f);
 reb_simulation_add_fmt(sim, "m a e inc Omega omega f", perturber_m, perturber_a, perturber_e, perturber_inc, perturber_Omega, perturber_omega, perturber_f);

 // Kill run if under Kozai inclination
 struct reb_orbit ob = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
 struct reb_orbit oc = reb_orbit_from_particle(sim->G, sim->particles[2], sim->particles[0]);
 angle = acos(reb_vec3d_dot(ob.hvec, oc.hvec) / (sqrt(reb_vec3d_length_squared(ob.hvec)) * sqrt(reb_vec3d_length_squared(oc.hvec))));
   
 if (angle * 180./M_PI < 40.){
   exit(1);
 }

  sim->integrator         = REB_INTEGRATOR_IAS15; // IAS15 is used for its adaptive timestep:
  sim->ri_ias15.adaptive_mode = 2;                // in a ZLK cycle the planet experiences close encounters during the high-eccentricity epochs.
  sim->heartbeat = heartbeat;                     // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation

  // Add REBOUNDx effects
  struct rebx_extras* rebx = rebx_attach(sim);
  struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
  rebx_add_force(rebx, effect);


  // Star is roughly Sun-like
  const double solar_k2 = 0.03;
  rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
  rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.28 * star.m * star.r * star.r);

  const double solar_spin_period = 25. * 2. * M_PI / 365.;
  const double solar_spin = (2 * M_PI) / solar_spin_period;
  rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default
  rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", 1e-8);

  // Planet
  planet_k2 = 0.5;
  rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
  rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * planet_m * planet_r * planet_r);
  struct reb_vec3d Omega_sv = reb_vec3d_mul(reb_vec3d_normalize(ob.hvec), pseudo_sync(planet_e, ob.n));
  rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);
  rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1e-5*(2*M_PI));


  // add GR precession:
  struct rebx_force* gr = rebx_load_force(rebx, "gr_potential");
  rebx_add_force(rebx, gr);
  rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units
  reb_simulation_move_to_com(sim);
  rebx_spin_initialize_ode(rebx, effect);

  // Set initial planet radius
  double lum = Etide(sim, rebx, &sim->particles[1], &sim->particles[0]);
  double rad = interpolate_mpb(lum);
  sim->particles[1].r = rad;


  reb_simulation_integrate(sim, tmax);
  rebx_free(rebx);
  reb_simulation_free(sim);

}

void heartbeat(struct reb_simulation* sim){
   // Radius inflation
   if(reb_simulation_output_check(sim, tscale)){
     struct rebx_extras* const rebx = sim->extras;
     struct reb_particle* p1 = &sim->particles[1];
     double lum = Etide(sim, rebx, &sim->particles[1], &sim->particles[0]);
     double rad = interpolate_mpb(lum);
     p1->r = rad;

     struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
     }
   }
}
