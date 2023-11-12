
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax;
int ind;
int stable=1;

char title[100] = "output_";
char title_stats[100] = "output_stability_stats";
char title_remove[100] = "rm -v output_";

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    sim->integrator         = REB_INTEGRATOR_WHFAST;
    sim->heartbeat          = heartbeat;

    ind = 0;
    if (argc == 2){
      strcat(title, argv[1]);
      strcat(title_remove, argv[1]);
      ind = atoi(argv[1]);
    }

    sim->rand_seed = ind;

    // Initial conditions
    // Santerne et al 2019
    struct reb_particle star = {0};
    star.m  = reb_random_uniform(sim, 1.16 - 0.04, 1.16 + 0.04);
    star.r = reb_random_uniform(sim, 1.273 - 0.015, 1.273 + 0.015) * 0.00465;
    reb_add(sim, star);

    // Planets
    double mearth = 3e-6;
    double ri = 0.05 * M_PI/180.;

    // b
    double mb = reb_random_uniform(sim, 6.89 - 0.88, 6.89 + 0.88) * mearth;
    double eb = reb_random_uniform(sim, 0.07 - 0.06, 0.07 + 0.06);
    double ab = reb_random_uniform(sim, 0.1283 - 1.5e-3, 0.1283 + 1.5e-3);
    double ib = reb_random_rayleigh(sim, ri);
    double Mb = reb_random_uniform(sim, 0, 2 * M_PI);

    double mc = reb_random_uniform(sim, 4.4 - 1.1, 4.4 + 1.1) * mearth;
    double ec = reb_random_uniform(sim, 0.04 - 0.03, 0.04 + 0.04);
    double ac = reb_random_uniform(sim, 0.2061 - 2.4e-3, 0.2061 + 2.4e-3);
    double ic = reb_random_rayleigh(sim, ri);
    double Mc = reb_random_uniform(sim, 0, 2 * M_PI);

    double md = reb_random_uniform(sim, 1.0, 4.6) * mearth;
    double ed = reb_random_uniform(sim, 0.06 - 0.06, 0.06 + 0.06);
    double ad = reb_random_uniform(sim, 0.88 - 0.01, 0.88 + 0.01);
    double id = reb_random_rayleigh(sim, ri);
    double Md = reb_random_uniform(sim, 0, 2 * M_PI);

    double me = reb_random_uniform(sim, 12. - 5., 12. + 5.) * mearth;
    double ee = reb_random_uniform(sim, 0.14 - 0.09, 0.14 + 0.09);
    double ae = reb_random_uniform(sim, 1.06 - 0.02, 1.06 + 0.03);
    double ie = reb_random_rayleigh(sim, ri);
    double Me = reb_random_uniform(sim, 0, 2 * M_PI);

    // This is the one we care abotu
    double mf = reb_random_uniform(sim, 12. - 3., 12. + 3.) * mearth;
    double ef = reb_random_uniform(sim, 0.004 - 0.003, 0.004 +0.009);
    double af = reb_random_uniform(sim, 1.37 - 0.02, 1.37 + 0.02);
    double incf = reb_random_rayleigh(sim, ri);
    double Mf = reb_random_uniform(sim, 0, 2 * M_PI);

    double planet_as[100] = {ab, ac, ad, ae, af};

    reb_add_fmt(sim, "m a e inc M", mb, ab, eb, ib, Mb);
    reb_add_fmt(sim, "m a e inc M", mc, ac, ec, ic, Mc);
    reb_add_fmt(sim, "m a e inc M", md, ad, ed, id, Md);
    reb_add_fmt(sim, "m a e inc M", me, ae, ee, ie, Me);
    reb_add_fmt(sim, "m a e inc M", mf, af, ef, incf, Mf);

    struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    sim->dt = o.P / 15.6789;

    // GR
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gr = rebx_load_force(rebx, "gr");
    rebx_add_force(rebx, gr);
    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32);

    reb_move_to_com(sim);

    struct reb_vec3d newz = reb_tools_angular_momentum(sim);
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    reb_simulation_irotate(rebx, rot);

    FILE* of = fopen(title, "w");
    fprintf(of, "t,nx,ny,nz\n");
    fclose(of);

    tmax = 1e5 * 2 * M_PI;
    reb_integrate(sim, tmax);

    for (unsigned i = 0; i < 5; i++){
      struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[i + 1], sim->particles[0]);
      if (fabs(orb.a - planet_as[i]) / planet_as[i] > 0.1){
        //printf("Unstable %d %f\n", i+1,orb.a);
        stable=0;
        system(title_remove);
        break;
      }
    }

    FILE* sf = fopen(title_stats, "a");
    fprintf(sf, "%d,%f,%f,%f,%f,%f,%d\n",mb/mearth,mc/mearth,md/mearth,me/mearth,mf/mearth,stable);
    fclose(sf);

    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 100. * 2 * M_PI)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen(title, "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      struct reb_orbit op = reb_tools_particle_to_orbit(sim->G, sim->particles[5], sim->particles[0]); // planet orbit
      struct reb_vec3d n = op.hvec;
      fprintf(of, "%f,%e\n", sim->t, n.x);
      fclose(of);
    }

    //if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
    //    reb_output_timing(sim, tmax);
    //}
}
