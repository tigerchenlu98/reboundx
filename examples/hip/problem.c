
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
int first_set=1;
int ntest=3;
double mf;
double me;

double planet_as[10] = {0.1283,0.2061,0.88,1.06,1.37};
double planet_aerrs[10] = {1.5e-3, 2.4e-3, 0.01, 0.03, 0.02};

char title[100] = "127mt_";
char title_stats[100] = "mt_stats";
char title_remove[100] = "rm -v 127mt_";

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();
    sim->integrator         = REB_INTEGRATOR_WHFAST;
    sim->heartbeat          = heartbeat;

    ind = 0;
    //double ob = 0;
    if (argc == 2){
      strcat(title, argv[1]);
      strcat(title_remove, argv[1]);
      ind = atoi(argv[1]);
      //ob = atoi(argv[2]);
    }

    sim->rand_seed = ind;
    double delta = reb_random_uniform(sim, 1.05, 1.15);

    // Initial conditions
    // Santerne et al 2019
    struct reb_particle star = {0};
    star.m  = 1.16;
    star.r = 1.273 * 0.00465;
    reb_simulation_add(sim, star);

    // Planets
    double mearth = 3e-6;
    double ri = 1.0 * M_PI/180.;

    // b
    double mb = 6.89 * mearth;
    double eb = 0.07;
    double ab = 0.1283;
    double ib = reb_random_rayleigh(sim, ri);
    double tb = reb_random_uniform(sim, 0, 2 * M_PI);
    //double Mb = reb_random_uniform(sim, 0, 2 * M_PI);

    double mc = 4.4 * mearth;
    double ec = 0.04;
    double ac = 0.2061;
    double ic = reb_random_rayleigh(sim, ri);
    double tc = reb_random_uniform(sim, 0, 2 * M_PI);
    //double Mc = reb_random_uniform(sim, 0, 2 * M_PI);

    double md = 4.6 * mearth;
    double ed = 0.06;
    double ad = 0.88 * reb_random_uniform(sim, 2., 3.);
    double id = reb_random_rayleigh(sim, ri);
    double td = reb_random_uniform(sim, 0, 2 * M_PI);
    //double Md = reb_random_uniform(sim, 0, 2 * M_PI);

    me = reb_random_uniform(sim, 12. - 5., 9.) * mearth;
    double ee = 0.14;
    double ae = ad * pow(4./3., 2./3.) * delta;//1.06;//reb_random_uniform(sim, 1.06 - 0.02, 1.06 + 0.03);
    //printf("%f %f\n", ad, ae);
    //exit(1);
    double ie = reb_random_rayleigh(sim, ri);
    double te = reb_random_uniform(sim, 0, 2 * M_PI);
    //double Me = reb_random_uniform(sim, 0, 2 * M_PI);

    // This is the one we care abotu
    mf = reb_random_uniform(sim, 9., 12.5) * mearth;
    double rho = 1.0 * pow(1.496e13, 3.) / (1.989e33); // 1 g/cm3 to rebound units
    double rf = pow(((3. * mf) / (4. * M_PI * rho)), 1./3.);
    double ef = 0.004;
    double af = ae * pow(3./2.,2./3.) * delta;//1.37;//reb_random_uniform(sim, 1.37 - 0.02, 1.37 + 0.02);
    double incf = reb_random_rayleigh(sim, ri);
    double tf = reb_random_uniform(sim, 0, 2 * M_PI);
    //double Mf = reb_random_uniform(sim, 0, 2 * M_PI);

    //double rhof = 1.0 * pow(1.496e13,3) / (1.989e33); // 1 g/cm^3 to solar masses/AU^3
    //rf =  pow((3 * mf / (4 * M_PI * rhof)), 1./3.);

    double planet_as[10] = {ab, ac, ad, ae, af};

    reb_simulation_add_fmt(sim, "primary m a e inc theta", star, mb, ab, eb, ib, tb);
    reb_simulation_add_fmt(sim, "primary m a e inc theta", star, mc, ac, ec, ic, tc);
    reb_simulation_add_fmt(sim, "primary m a e inc theta", star, md, ad, ed, id, td);
    reb_simulation_add_fmt(sim, "primary m a e inc theta", star, me, ae, ee, ie, te);
    reb_simulation_add_fmt(sim, "primary m r a e inc theta", star, mf, rf, af, ef, incf, tf);

    //double rin = 1.5 * rf;
    //double rho_ring = 0.5;
    //double rout = 2.45 * pow(1./rho_ring, 1./3.) * rf;
    //for (unsigned int i = 0; i < ntest; i++){
    //reb_add_fmt(sim, "primary a inc Omega", sim->particles[2], rout, 30.*M_PI/180.);
    //}

    // tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* mof = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mof);

    double tau_0 = -1e7 * 2 * M_PI;
    double kappa = 100.;
    double beta = -1.7;
    for (unsigned int i = 3; i < sim->N; i++){
        double tau_ai = tau_0 * pow(planet_as[i-1], beta);
        double tau_ei = tau_ai / kappa;
        rebx_set_param_double(rebx, &sim->particles[i].ap, "tau_a", tau_ai);
        rebx_set_param_double(rebx, &sim->particles[i].ap, "tau_e", tau_ei);
    }

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    double planet_k2 = 0.6;
    rebx_set_param_double(rebx, &sim->particles[5].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[5].ap, "I", 0.25 * mf * rf * rf);

    const double spin_period_p = ((10. / 24.) / 365.) * 2. * M_PI; // hours to reb years
    const double spin_p = (2. * M_PI) / spin_period_p;
    const double theta_p = 1. * M_PI / 180.;
    const double phi_p = 180. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[5].ap, "Omega", Omega_sv);

    const double planet_Q = 1e5;
    struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[5], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[5].ap, "tau", 1./(2.*orb.n*planet_Q));


    struct reb_vec3d newz = reb_vec3d_add(reb_simulation_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    if (isnan(rot.r)) {
      rot = reb_rotation_identity();
    }
    reb_simulation_move_to_com(sim);
    rebx_simulation_irotate(rebx, rot);
    rebx_spin_initialize_ode(rebx, effect);

    //FILE* of = fopen(title, "w");
    //fprintf(of, "t,magp,thetap,phip,Omegap\n");
    //fclose(of);

    system(title_remove);
    //FILE* of = fopen(title, "w");
    //fprintf(of, "t,mag,theta,phi,sx,sy,sz,ad,ae,af,omega,mf\n");
    //fprintf(of, "t,inc,Omega,nx\n");
    //for (unsigned int i = 0; i < ntest; i++){
    //fprintf(of, ",at,it");
    //}
    //fprintf(of, "\n");
    //fclose(of);

    struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    tmax = 2e5*2*M_PI;//o.P * 1e8;
    sim->dt = o.P / 10.12345;
    reb_simulation_integrate(sim, tmax);

    //for (unsigned i = 0; i < 5; i++){
    struct reb_orbit orf = reb_orbit_from_particle(sim->G, sim->particles[5], sim->particles[0]);
    struct reb_vec3d n1 = orf.hvec;

    struct reb_particle* planet = &sim->particles[5];
    struct reb_vec3d* Omega_p_inv = rebx_get_param(rebx, planet->ap, "Omega");
    struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, n1);
    struct reb_rotation prot = reb_rotation_init_to_new_axes(n1, line_of_nodes); // Arguments to this function are the new z and x axes
    if (isnan(prot.r)) {
      prot = reb_rotation_identity();
    }
    struct reb_vec3d srot = reb_vec3d_rotate(*Omega_p_inv, prot);

    double magp;
    double thetap;
    double phip;
    reb_tools_xyz_to_spherical(srot, &magp, &thetap, &phip);

    FILE* sf = fopen(title_stats, "a");
    fprintf(sf, "%d %f\n", ind, thetap * 180./M_PI);
    fclose(sf);
        //stable = 0;
        //system(title_remove);
        //break;
      //}
    //}
/*
    if (stable){
      FILE* sf = fopen(title_stats, "a");
      fprintf(sf, "Stable %d %f %f %f %f %f\n", ind, ob.a, oc.a, od.a, oe.a, orf.a);
      fclose(sf);
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
/*
    if(reb_simulation_output_check(sim, 1. * 2 * M_PI)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* planet = &sim->particles[5];

      //struct reb_orbit ob = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
      //struct reb_orbit oc = reb_orbit_from_particle(sim->G, sim->particles[2], sim->particles[0]);
      struct reb_orbit od = reb_orbit_from_particle(sim->G, sim->particles[3], sim->particles[0]);
      struct reb_orbit oe = reb_orbit_from_particle(sim->G, sim->particles[4], sim->particles[0]);
      struct reb_orbit orf = reb_orbit_from_particle(sim->G, sim->particles[5], sim->particles[0]);
      struct reb_vec3d n1 = orf.hvec;

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
      fprintf(sf, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",sim->t,mag_p,theta_p,phi_p,srot.x,srot.y,srot.z,od.a,oe.a,orf.a,orf.Omega,sim->particles[5].m);
      //fprintf(sf, "%f,%f,%f,%f\n",sim->t,o1.inc,o1.Omega,n1.x);
      fclose(sf);


    }
    */


    //if(reb_simulation_output_check(sim, 10.)){        // outputs to the screen
    //    reb_simulation_output_timing(sim, tmax);
    //}

    struct reb_orbit orbd = reb_orbit_from_particle(sim->G, sim->particles[3], sim->particles[0]);
    if (orbd.a < 0.881 && first_set){
      for (unsigned int i = 3; i < sim->N; i++){
          rebx_set_param_double(sim->extras, &sim->particles[i].ap, "tau_a", INFINITY);
          rebx_set_param_double(sim->extras, &sim->particles[i].ap, "tau_e", INFINITY);
      }
      first_set = 0;
    }

}
