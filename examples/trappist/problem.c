#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "tides_spin.c"

char title[100] = "output_4.txt";
void heartbeat(struct reb_simulation* r);
double tmax = 5000. * M_PI * 2.;
int planets = 7;

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
   star.m  = 0.0898;
   star.r = 0.1192 * 0.00465;
   reb_add(sim, star);

   // Planets
   double m1 = 4.596e-5 * star.m;
   double r1 = 0.08590 * star.r;
   double a1 = 20.843 * star.r;
   double e1 = 0.003055;
   double inc1 = 89.728 * M_PI/180.;
   double pomega1 = 134.73 * M_PI/180.;
   double Omega1 = 0.;
   double f1 = 45.98212 * M_PI/180.;

   double m2 = 4.374e-5 * star.m;
   double r2 = 0.08440 * star.r;
   double a2 = 28.549 * star.r;
   double e2 = 0.000550;
   double inc2 = 89.778 * M_PI/180.;
   double pomega2 = 1.041627 * M_PI/180.;
   double Omega2 = 0.;
   double f2 = -8.568899 * M_PI/180.;

   double m3 = 1.297e-5 * star.m;
   double r3 = 0.06063 * star.r;
   double a3 = 40.216 * star.r;
   double e3 = 0.005633;
   double inc3 = 89.896 * M_PI/180.;
   double pomega3 = 151.706133229934 * M_PI/180.;
   double Omega3 = 0.;
   double f3 = 15.0620262596043 * M_PI/180.;

   double m4 = 2.313e-5 * star.m;
   double r4 = 0.07079 * star.r;
   double a4 = 52.855 * star.r;
   double e4 = 0.006325;
   double inc4 = 89.793 * M_PI/180.;
   double pomega4 = 313.206087730948 * M_PI/180.;
   double Omega4 = 0.;
   double f4 = -217.102135813412 * M_PI/180.;

   double m5 = 3.475e-5 * star.m;
   double r5 = 0.08040 * star.r;
   double a5 = 69.543 * star.r;
   double e5 = 0.008415;
   double inc5 = 89.740 * M_PI/180.;
   double pomega5 = 183.474407367532 * M_PI/180.;
   double Omega5 = 0.;
   double f5 = -59.9711817502667 * M_PI/180.;

   double m6 = 4.418e-5 * star.m;
   double r6 = 0.08692 * star.r;
   double a6 = 84.591 * star.r;
   double e6 = 0.004010;
   double inc6 = 89.742 * M_PI/180.;
   double pomega6 = 18.615692007331 * M_PI/180.;
   double Omega6 = 0.;
   double f6 = 77.6950168190991 * M_PI/180.;

   double m7 = 4.596e-5 * star.m;
   double r7 = 0.08590 * star.r;
   double a7 = 20.843 * star.r;
   double e7 = 0.003055;
   double inc7 = 89.728 * M_PI/180.;
   double pomega7 = 180.313946334794 * M_PI/180.;
   double Omega7 = 0.;
   double f7 = 69.3197659604004 * M_PI/180.;

   double ms[7] = {m1,m2,m3,m4,m5,m6,m7};
   double rs[7] = {r1,r2,r3,r4,r5,r6,r7};

   reb_add_fmt(sim, "m r a e inc pomega Omega f", m1, r1, a1, e1, inc1, pomega1, Omega1, f1);
   reb_add_fmt(sim, "m r a e inc pomega Omega f", m2, r2, a2, e2, inc2, pomega2, Omega2, f2);
   reb_add_fmt(sim, "m r a e inc pomega Omega f", m3, r3, a3, e3, inc3, pomega3, Omega3, f3);
   reb_add_fmt(sim, "m r a e inc pomega Omega f", m4, r4, a4, e4, inc4, pomega4, Omega4, f4);
   reb_add_fmt(sim, "m r a e inc pomega Omega f", m5, r5, a5, e5, inc5, pomega5, Omega5, f5);
   reb_add_fmt(sim, "m r a e inc pomega Omega f", m6, r6, a6, e6, inc6, pomega6, Omega6, f6);
   reb_add_fmt(sim, "m r a e inc pomega Omega f", m7, r7, a7, e7, inc7, pomega7, Omega7, f7);

   struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
   sim->dt = orb.P / 15.6789;     // initial timestep as a function of orbital period

   // Add REBOUNDx effects
   // First tides_spin
   struct rebx_extras* rebx = rebx_attach(sim);
   struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
   rebx_add_force(rebx, effect);

   // Sun - fiducial values, should not impact simulation significantly
   const double solar_k2 = 0.07;
   rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
   rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * star.m * star.r * star.r);
   const double solar_spin_period = 27. * 2. * M_PI / 365.;
   const double solar_spin = (2 * M_PI) / solar_spin_period;
   rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.y=solar_spin}); // Omega_x = Omega_y = 0 by default
   rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", 1./(2.*1000000.*orb.n));

   // Planets
   double k2 = 0.299;
   double tau_p = 712.37 / (86400. * 365.25 * 2 * M_PI); // Bolmont et al 2015, seconds to reb years
   double spin_rate = 0.2617 * 24. * 365.25 * 2 * M_PI; // rad/hr to rad/reb years
   for (unsigned int i = 0; i < planets; i++){
     rebx_set_param_double(rebx, &sim->particles[i+1].ap, "k2", k2);
     rebx_set_param_double(rebx, &sim->particles[i+1].ap, "I", 0.3 * ms[i] * rs[i] * rs[i]);

     struct reb_orbit op = reb_tools_particle_to_orbit(sim->G, sim->particles[i+1], sim->particles[0]);
     struct reb_vec3d spin_vec = {};
     spin_vec.x = op.hvec.x * spin_rate;
     spin_vec.y = op.hvec.y * spin_rate;
     spin_vec.z = op.hvec.z * spin_rate;

     rebx_set_param_vec3d(rebx, &sim->particles[i+1].ap, "Omega", spin_vec);
     rebx_set_param_double(rebx, &sim->particles[i+1].ap, "tau", tau_p);
   }

   // Set obliquity of the 4th planet
   double theta_e = 23.43389 * M_PI/180.;
   struct reb_vec3d spin_vec_e = reb_tools_spherical_to_xyz(spin_rate, inc5 - theta_e, 90.0 * M_PI/180.);
   rebx_set_param_vec3d(rebx, &sim->particles[4].ap, "Omega", spin_vec_e);

   // General Relativity
   //struct rebx_force* gr = rebx_load_force(rebx, "gr");
   //rebx_add_force(rebx, gr);
   //rebx_set_param_double(rebx, &gr->ap, "c", 10065.32);

   // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
   // nodes with the invariable plane (along z cross newz)
   struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
   struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
   struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);

   rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
   reb_move_to_com(sim);
   rebx_spin_initialize_ode(rebx, effect);


   reb_integrate(sim, tmax);
   rebx_free(rebx);
   reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
   // Output spin and orbital information to file
   if(reb_output_check(sim, 0.001 / (365.25 * 2 * M_PI))){        // outputs every 10 REBOUND years
     struct rebx_extras* const rebx = sim->extras;
     FILE* of = fopen(title, "a");
     if (of==NULL){
         reb_error(sim, "Can not open file.");
         return;
     }

     // Planet e spin vector
     struct reb_particle* e = &sim->particles[4];
     struct reb_vec3d* Omega_e = rebx_get_param(rebx, e->ap, "Omega");
     fprintf(of, "%f,%f,%f,%f",sim->t, Omega_e->x, Omega_e->y, Omega_e->z);

     for (unsigned int i = 0; i < planets; i++){
       struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[i+1], sim->particles[0]);
       fprintf(of, ",%f,%f,%f,%f,%f,%f", o.a, o.e, o.inc, o.Omega, o.omega, o.f); // print spins and orbits
     }
     fprintf(of,"\n");
     fclose(of);
   }
}
