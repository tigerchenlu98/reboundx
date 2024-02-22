
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
double mb;
double mc;
double md;
double cJ2;
double dJ2;

//char title[100] = "ls";
char title_stats[100] = "222_kepler_mt_stats";
//char title_remove[100] = "rm -v ls";

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();
    sim->integrator         = REB_INTEGRATOR_WHFAST;
    sim->heartbeat          = heartbeat;

    if (argc == 2){
      ind = atoi(argv[1]);
    }

    sim->rand_seed = ind;
    double delta = 1.02;


    // Initial conditions
    // Santerne et al 2019
    struct reb_particle star = {0};
    star.m  = 0.985;
    reb_simulation_add(sim, star);

    // Planets
    double mearth = 3e-6;
    double ri = 0.5 * M_PI/180.;

    // b
    mb = reb_random_uniform(sim, 2.1-0.6, 2.1+1.5) * mearth;
    double eb = 0.0;
    double ab = 0.2514;
    double ib = reb_random_rayleigh(sim, ri);
    double tb = reb_random_uniform(sim, 0, 2 * M_PI);

    mc = reb_random_uniform(sim, 4.0-0.4, 4.0+0.4) * mearth;
    double rho_c = reb_random_uniform(sim, 1.0, 5.5);
    double rho_c_converted = rho_c * pow(1.496e13, 3.) / (1.989e33); // 5 g/cm3 to rebound units
    double rc = pow(((3. * mc) / (4. * M_PI * rho_c_converted)), 1./3.);
    double ec = 0.0;
    double ac = 0.384 * 3.0;
    double ic = reb_random_rayleigh(sim, ri);
    double tc = reb_random_uniform(sim, 0, 2 * M_PI);

    md = reb_random_uniform(sim, 7.6-1.1, 7.6+1.1) * mearth;
    double rho_d = reb_random_uniform(sim, 1.0, 5.5);
    double rho_d_converted = rho_d * pow(1.496e13, 3.) / (1.989e33); // 5 g/cm3 to rebound units
    double rd = pow(((3. * md) / (4. * M_PI * rho_d_converted)), 1./3.);
    double ed = 0.0;
    double ad = ac * pow(3./2., 2./3.) * delta;
    double id = reb_random_rayleigh(sim, ri);
    double td = reb_random_uniform(sim, 0, 2 * M_PI);

    double planet_as[10] = {ab, ac, ad};

    reb_simulation_add_fmt(sim, "primary m a e inc theta", star, mb, ab, eb, ib, tb);
    reb_simulation_add_fmt(sim, "primary m r a e inc theta", star, mc, rc, ac, ec, ic, tc);
    reb_simulation_add_fmt(sim, "primary m r a e inc theta", star, md, rd, ad, ed, id, td);

    // tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* mof = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mof);

    double tau_0 = -5e7 * 2 * M_PI;
    double kappa = 100.;
    double beta = -1.7;
    for (unsigned int i = 2; i < sim->N; i++){
        double tau_ai = tau_0 * pow(planet_as[i-1], beta);
        double tau_ei = tau_ai / kappa;
        rebx_set_param_double(rebx, &sim->particles[i].ap, "tau_a", tau_ai);
        rebx_set_param_double(rebx, &sim->particles[i].ap, "tau_e", tau_ei);
    }

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    const double planet_Q = 1e5;
    const double theta_p = 1. * M_PI / 180.;
    const double phi_p = 180. * M_PI / 180;

    double ck2 = reb_random_uniform(sim, 0.1, 0.8);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", ck2);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "I", 0.25 * mc * rc * rc);

    const double c_spin_hrs = reb_random_uniform(sim, 5./24., 1.);
    const double c_spin_period = (c_spin_hrs / 365.) * 2. * M_PI; // hours to reb years
    const double c_spin_p = (2. * M_PI) / c_spin_period;
    struct reb_vec3d c_Omega = reb_tools_spherical_to_xyz(c_spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[2].ap, "Omega", c_Omega);

    cJ2 = (c_spin_p * c_spin_p * rc * rc * rc) / (3. * sim->G * mc) * ck2;
    double corbn = sqrt(sim->G * (star.m + mc) / (0.384 * 0.384 * 0.384));
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau", 1./(2.*corbn*planet_Q));

    double dk2 = reb_random_uniform(sim, 0.1, 0.8);
    rebx_set_param_double(rebx, &sim->particles[3].ap, "k2", dk2);
    rebx_set_param_double(rebx, &sim->particles[3].ap, "I", 0.25 * md * rd * rd);

    const double d_spin_hrs = reb_random_uniform(sim, 5./24., 1.);
    const double d_spin_period = (d_spin_hrs / 365.) * 2. * M_PI; // hours to reb years
    const double d_spin_p = (2. * M_PI) / d_spin_period;
    struct reb_vec3d d_Omega = reb_tools_spherical_to_xyz(d_spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[3].ap, "Omega", d_Omega);

    dJ2 = (d_spin_p * d_spin_p * rd * rd * rd) / (3. * sim->G * md) * dk2;
    double dorbn = sqrt(sim->G * (star.m + md) / (0.509 * 0.509 * 0.509));
    rebx_set_param_double(rebx, &sim->particles[3].ap, "tau", 1./(2.*dorbn*planet_Q));


    struct reb_vec3d newz = reb_vec3d_add(reb_simulation_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    //struct reb_vec3d newz = reb_simulation_angular_momentum(sim);
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    if (isnan(rot.r)) {
      rot = reb_rotation_identity();
    }

    reb_simulation_move_to_com(sim);
    rebx_simulation_irotate(rebx, rot);
    rebx_spin_initialize_ode(rebx, effect);

    const double c_alpha_init = 0.5 * (star.m / mc) * pow((rc / ac), 3.) * (ck2 / 0.25) * c_spin_p;
    const double d_alpha_init = 0.5 * (star.m / md) * pow((rd / ad), 3.) * (dk2 / 0.25) * d_spin_p;

    struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    tmax = 1e3*2*M_PI;//o.P * 1e8;
    sim->dt = o.P / 15.12345;
    reb_simulation_integrate(sim, tmax);

    // Spin axes
    struct reb_orbit orbb = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    struct reb_orbit orbc = reb_orbit_from_particle(sim->G, sim->particles[2], sim->particles[0]);
    struct reb_orbit orbd = reb_orbit_from_particle(sim->G, sim->particles[3], sim->particles[0]);
    struct reb_vec3d nc = orbc.hvec;
    struct reb_vec3d nd = orbd.hvec;

    struct reb_particle* pc = &sim->particles[2];
    struct reb_vec3d* cOmega_inv = rebx_get_param(rebx, pc->ap, "Omega");
    struct reb_vec3d cline_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, nc);
    struct reb_rotation crot = reb_rotation_init_to_new_axes(nc, cline_of_nodes); // Arguments to this function are the new z and x axes
    if (isnan(crot.r)) {
      crot = reb_rotation_identity();
    }
    struct reb_vec3d cOmega = reb_vec3d_rotate(*cOmega_inv, crot);

    double cmagp;
    double cthetap;
    double cphip;
    reb_tools_xyz_to_spherical(cOmega, &cmagp, &cthetap, &cphip);

    struct reb_particle* pd = &sim->particles[3];
    struct reb_vec3d* dOmega_inv = rebx_get_param(rebx, pd->ap, "Omega");
    struct reb_vec3d dline_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, nd);
    struct reb_rotation drot = reb_rotation_init_to_new_axes(nd, dline_of_nodes); // Arguments to this function are the new z and x axes
    if (isnan(drot.r)) {
      drot = reb_rotation_identity();
    }
    struct reb_vec3d dOmega = reb_vec3d_rotate(*dOmega_inv, drot);

    double dmagp;
    double dthetap;
    double dphip;
    reb_tools_xyz_to_spherical(dOmega, &dmagp, &dthetap, &dphip);

    const double c_alpha_final = 0.5 * (star.m/mc) * pow((rc / orbc.a), 3.) * (ck2 / 0.25) * cmagp;
    const double d_alpha_final = 0.5 * (star.m/md) * pow((rd / orbd.a), 3.) * (dk2 / 0.25) * dmagp;

    const double new_cJ2 = (cmagp * cmagp * cmagp * rc * rc * rc) / (3. * sim->G * mc) * ck2;
    const double cLR = pow(new_cJ2 * rc * rc * orbc.a * orbc.a * orbc.a * pow((1. - orbc.e * orbc.e),3./2.) * (mc/star.m), 1./5.);

    const double new_dJ2 = (dmagp * dmagp * dmagp * rd * rd * rd) / (3. * sim->G * md) * dk2;
    const double dLR = pow(new_dJ2 * rd * rd * orbd.a * orbd.a * orbd.a * pow((1. - orbd.e * orbd.e),3./2.) * (md/star.m), 1./5.);

    FILE* sf = fopen(title_stats, "a");
    fprintf(sf, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%e,%e,%e,%e,%f,%f,%f\n", ind, rho_c, rho_d, mb/mearth, mc/mearth,md/mearth, rc/4.259e-5, rd/4.259e-5, cthetap * 180./M_PI, dthetap * 180./M_PI,ck2,dk2,c_spin_hrs * 24.,d_spin_hrs * 24.,cLR/rc,dLR/rd,cJ2,dJ2,c_alpha_init,d_alpha_init,orbb.a,orbc.a,orbd.a);
    fclose(sf);
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
      //struct reb_orbit od = reb_orbit_from_particle(sim->G, sim->particles[3], sim->particles[0]);
      //struct reb_orbit oe = reb_orbit_from_particle(sim->G, sim->particles[4], sim->particles[0]);
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
      //fprintf(sf, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",sim->t,mag_p,theta_p,phi_p,srot.x,srot.y,srot.z,od.a,oe.a,orf.a,orf.Omega,sim->particles[5].m);
      fprintf(sf, "%f,%f,%f,%f\n",sim->t,orf.inc,orf.Omega,n1.x);
      fclose(sf);


    }



    if(reb_simulation_output_check(sim, 10.)){        // outputs to the screen
        reb_simulation_output_timing(sim, tmax);
    }
*/
    if (first_set){
      struct reb_orbit orbc = reb_orbit_from_particle(sim->G, sim->particles[3], sim->particles[0]);
      if (orbc.a < 0.509){
        for (unsigned int i = 2; i < sim->N; i++){
            rebx_set_param_double(sim->extras, &sim->particles[i].ap, "tau_a", INFINITY);
            rebx_set_param_double(sim->extras, &sim->particles[i].ap, "tau_e", INFINITY);
        }
        first_set = 0;
      }
    }

}
