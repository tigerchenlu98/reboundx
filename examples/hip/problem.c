
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
int ntest=3;
double mf;
double me;

double planet_as[10] = {0.1283,0.2061,0.88,1.06,1.37};
double planet_aerrs[10] = {1.5e-3, 2.4e-3, 0.01, 0.03, 0.02};

char title[100] = "ob_evolution_";
//char title_stats[100] = "stability_stats";
//char title_remove[100] = "rm -v low_ob_";

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();
    sim->integrator         = REB_INTEGRATOR_WHFAST;
    sim->heartbeat          = heartbeat;

    ind = 0;
    if (argc == 2){
      strcat(title, argv[1]);
      //strcat(title_remove, argv[1]);
      ind = atoi(argv[1]);
    }

    sim->rand_seed = ind;

    // Initial conditions
    // Santerne et al 2019
    struct reb_particle star = {0};
    star.m  = 1.16;
    star.r = 1.273 * 0.00465;
    reb_simulation_add(sim, star);

    // Planets
    double mearth = 3e-6;
    double ri = 0.05 * M_PI/180.;

    // b
    double mb = 6.89 * mearth;
    double eb = 0.07;
    double ab = 0.1283;
    double ib = reb_random_uniform(sim, -1. * ri, ri);
    double Mb = reb_random_uniform(sim, 0, 2 * M_PI);

    double mc = 4.4 * mearth;
    double ec = 0.04;
    double ac = 0.2061;
    double ic = reb_random_uniform(sim, -1. * ri, ri);
    double Mc = reb_random_uniform(sim, 0, 2 * M_PI);

    double md = 4.6 * mearth;
    double ed = 0.06;
    double ad = 0.88;
    double id = 0.0;//reb_random_uniform(sim, -1. * ri, ri);
    double Md = reb_random_uniform(sim, 0, 2 * M_PI);

    me = reb_random_uniform(sim, 12. - 5., 12. + 5.) * mearth;
    double ee = 0.14;
    double ae = 1.06;//reb_random_uniform(sim, 1.06 - 0.02, 1.06 + 0.03);
    double ie = reb_random_uniform(sim, -1. * ri, ri);
    double Me = reb_random_uniform(sim, 0, 2 * M_PI);

    // This is the one we care abotu
    mf = reb_random_uniform(sim, 12. - 3., 12. + 3.) * mearth;
    double rho = 1.0 * pow(1.496e13, 3.) / (1.989e33); // 1 g/cm3 to rebound units
    double rf = pow(((3. * mf) / (4. * M_PI * rho)), 1./3.);
    double ef = 0.004;
    double af = 1.37;//reb_random_uniform(sim, 1.37 - 0.02, 1.37 + 0.02);
    double incf = reb_random_uniform(sim, -1. * ri, ri);
    double Mf = reb_random_uniform(sim, 0, 2 * M_PI);

    //double rhof = 1.0 * pow(1.496e13,3) / (1.989e33); // 1 g/cm^3 to solar masses/AU^3
    //rf =  pow((3 * mf / (4 * M_PI * rhof)), 1./3.);

    double planet_as[10] = {ab, ac, ad, ae, af};

    reb_simulation_add_fmt(sim, "primary m a e inc M", star, mb, ab, eb, ib, Mb);
    reb_simulation_add_fmt(sim, "primary m a e inc M", star, mc, ac, ec, ic, Mc);
    //reb_simulation_add_fmt(sim, "m a e inc M", md, ad, ed, id, Md);
    reb_simulation_add_fmt(sim, "m a e inc M", me, ae, ee, ie, Me);
    reb_simulation_add_fmt(sim, "m r a e inc M", mf, rf, af, ef, incf, Mf);

    //double rin = 1.5 * rf;
    //double rho_ring = 0.5;
    //double rout = 2.45 * pow(1./rho_ring, 1./3.) * rf;
    //for (unsigned int i = 0; i < ntest; i++){
    //reb_add_fmt(sim, "primary a inc Omega", sim->particles[2], rout, 30.*M_PI/180.);
    //}

    // tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    double planet_k2 = 0.6;
    rebx_set_param_double(rebx, &sim->particles[4].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[4].ap, "I", 0.25 * mf * rf * rf);

    const double spin_period_p = (10. / 365.) * 2. * M_PI; // days to reb years
    const double spin_p = (2. * M_PI) / spin_period_p;
    const double theta_p = 5. * M_PI / 180.;
    const double phi_p = 0. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[4].ap, "Omega", Omega_sv);

    const double planet_Q = 1e5;
    struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[4], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[4].ap, "tau", 1./(2.*orb.n*planet_Q));

    reb_simulation_move_to_com(sim);

    struct reb_vec3d newz = reb_vec3d_add(reb_simulation_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    if (isnan(rot.r)) {
      rot = reb_rotation_identity();
    }
    rebx_simulation_irotate(rebx, rot);
    rebx_spin_initialize_ode(rebx, effect);

    FILE* of = fopen(title, "w");
    fprintf(of, "t,magp,thetap,phip\n");
    fclose(of);
/*
    system(title_remove);
    FILE* of = fopen(title, "w");
    fprintf(of, "t");
    //for (unsigned int i = 0; i < ntest; i++){
    fprintf(of, ",at,it");
    //}
    fprintf(of, "\n");
    fclose(of);
    */
    struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    tmax = 1e4*2*M_PI;//o.P * 1e8;
    sim->dt = o.P / 15.12345;
    reb_simulation_integrate(sim, tmax);
/*
    for (unsigned i = 0; i < 5; i++){
      struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[i + 1], sim->particles[0]);
      if (fabs(orb.a - planet_as[i]) / planet_as[i] > 0.1){
        //printf("Unstable %d %f\n", i+1,orb.a);
        stable=0;
        system(title_remove);
        break;
      }
    }
*/
/*
    FILE* sf = fopen(title_stats, "a");
    fprintf(sf, "%d,%f,%f,%d\n",ind,me,mf,1);
    fclose(sf);
*/
    rebx_free(rebx);
    reb_simulation_free(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file

    if(reb_simulation_output_check(sim, 1. * 2 * M_PI)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* planet = &sim->particles[4];

      struct reb_orbit o1 = reb_orbit_from_particle(sim->G, *planet, *sun);
      struct reb_vec3d n1 = o1.hvec;

      struct reb_vec3d* Omega_p_inv = rebx_get_param(rebx, planet->ap, "Omega");
      struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, n1);
      struct reb_rotation rot = reb_rotation_init_to_new_axes(n1, line_of_nodes); // Arguments to this function are the new z and x axes
      if (isnan(rot.r)) {
        rot = reb_rotation_identity();
      }
      struct reb_vec3d srot = reb_vec3d_rotate(*Omega_p_inv, rot);

      double mag_p;
      double theta_p;
      double phi_p;
      reb_tools_xyz_to_spherical(srot, &mag_p, &theta_p, &phi_p);

      FILE* sf = fopen(title, "a");
      fprintf(sf, "%f,%f,%f,%f\n",sim->t,mag_p,theta_p,phi_p);
      fclose(sf);


    }


    if(reb_simulation_output_check(sim, 10.)){        // outputs to the screen
        reb_simulation_output_timing(sim, tmax);
    }
}
