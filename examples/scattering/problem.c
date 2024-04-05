#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "rebound.h"
#include "particle.h"

void heartbeat(struct reb_simulation* r);

double MASSC = 3.06 * 9.55e-4;
double DELTA = 3.;
int ind;
int e_start;
double tmax = 1e7*2*M_PI;
int nbodies = 3;

double a1(double d){
  return (4.192 * (3.51376e17 + 1.09428e15 * d*d)) / ((3.42236e8 + 3.30799e7 * d) * (3.42236e8 + 3.30799e7 * d));
}

double a2(double d){
  return (-4.192 * (3.51376e17 + 1.09428e15 * d*d)) / (-1.17125e17 + 1.09428e15 * d * d);
}

double a3(double d){
  return (4.192 * (3.51376e17 + 1.09428e15 * d*d)) / ((-3.42236e8 + 3.30799e7 * d) * (-3.42236e8 + 3.30799e7 * d));
}

char title_stats[100] = "scattering_Delta3_stats";

int main(int argc, char* argv[]){
    // Initialize masses

    struct reb_simulation* r = reb_simulation_create();
    struct reb_particle star = {0};
    star.m = 1;
    reb_simulation_add(r, star);

    ind = 0;
    if (argc == 2){
      //strcat(title, argv[1]);
      //strcat(title_remove, argv[1]);
      ind = atoi(argv[1]);
    }

    r->rand_seed = ind;

    double smas[3] = {a1(DELTA), a2(DELTA), a3(DELTA)};

    for (int i = 0; i < nbodies; i++){
      reb_simulation_add_fmt(r, "m a e inc omega Omega f hash", MASSC, smas[i], reb_random_uniform(r, 0.01, 0.05), reb_random_uniform(r, 0.0, 2.*M_PI/180.), reb_random_uniform(r, 0.0, 2.*M_PI), reb_random_uniform(r, 0.0, 2.*M_PI), reb_random_uniform(r, 0.0, 2.*M_PI), i+1);
    }

    reb_simulation_move_to_com(r);

    r->integrator = REB_INTEGRATOR_IAS15;
    r->exit_max_distance = 1e5;
    r->track_energy_offset = 1;

    double e_init = reb_simulation_energy(r);

    while (r->t < tmax){
       int retval = reb_simulation_integrate(r, tmax);
       if (retval == REB_STATUS_ESCAPE){
          double Ei = reb_simulation_energy(r);

          // Find and remove the particle
          int remove_ind;
          for (int i = 1; i < nbodies+1; i++){

              struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
              if (p != NULL){
                double dx = p->x;
                double dy = p->y;
                double dz = p->z;
                double d2 = dx*dx+dy*dy+dz*dz;

                if (d2>r->exit_max_distance * r->exit_max_distance){
                    remove_ind = i;
                }
              }
          }
          reb_simulation_remove_particle_by_hash(r, remove_ind, 1);
          reb_simulation_move_to_com(r);

          r->energy_offset += Ei - reb_simulation_energy(r);

          // Only one planet left
          if (r->N == 2){
            FILE* tf = fopen(title_stats, "a");
            fprintf(tf, "%d,%d,%e", ind, r->N-1, fabs((reb_simulation_energy(r) - e_init)/e_init));
            struct reb_particle* sun = &r->particles[0];
            for (int i = 1; i < nbodies+1;i++){
              struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
              if (p != NULL){
                struct reb_orbit o = reb_orbit_from_particle(r->G, *p, *sun);
                fprintf(tf, ",%f,%f,%f", o.a, o.e, o.inc);
              }
              else{
                fprintf(tf, ",-1,-1,-1");
              }
            }
            fprintf(tf, "%e\n", r->t);
            fclose(tf);
            exit(1);
          }
       }
    }

    FILE* tf = fopen(title_stats, "a");
    fprintf(tf, "%d,%d,%e", ind, r->N-1, fabs((reb_simulation_energy(r) - e_init)/e_init));
    struct reb_particle* sun = &r->particles[0];
    for (int i = 1; i < nbodies+1;i++){
      struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
      if (p != NULL){
        struct reb_orbit o = reb_orbit_from_particle(r->G, *p, *sun);
        fprintf(tf, ",%f,%f,%f", o.a, o.e, o.inc);
      }
      else{
        fprintf(tf, ",-1,-1,-1");
      }
    }
    fprintf(tf, "%e\n", tmax);
    fclose(tf);

    reb_simulation_free(r);
}
