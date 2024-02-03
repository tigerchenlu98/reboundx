
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
double mf;
int ind;
double LR;
int first_check=1;

struct rebx_interpolator* sxfunc;
struct rebx_interpolator* syfunc;
struct rebx_interpolator* szfunc;
struct rebx_extras* rebx;

char title[100] = "spinp_test_";
// char title_stats[100] = "migration_stats";
char title_remove[100] = "rm -v spinp_test_";

double laplace_radius(struct reb_simulation* sim, struct rebx_extras* const rebx, struct reb_particle* planet, struct reb_particle* star){
  const double* k2 = rebx_get_param(rebx, planet->ap, "k2");
  const struct reb_vec3d* Omega = rebx_get_param(rebx, planet->ap, "Omega");

  double spin_p = sqrt(reb_vec3d_length_squared(*Omega));

  struct reb_orbit o = reb_orbit_from_particle(sim->G, *planet, *star);
  const double j2 = ((*k2) / 3.) * (spin_p * spin_p * planet->r * planet->r * planet->r) / (sim->G * planet->m);
  const double lr =  pow(2. * j2 * planet->r * planet->r * o.a * o.a * o.a * pow((1 - o.e * o.e),(3./2.)) * planet->m / star->m, (1./5.));
  return lr;
}

double laplace_equilibrium(struct reb_simulation* sim, struct rebx_extras* const rebx, struct reb_particle* planet, struct reb_particle* star, double d){
  struct reb_orbit o = reb_orbit_from_particle(sim->G, *planet, *star);
  struct reb_vec3d nvec = o.hvec;
  const struct reb_vec3d* Omega = rebx_get_param(rebx, planet->ap, "Omega");

  struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, nvec);
  struct reb_rotation rot = reb_rotation_init_to_new_axes(nvec, line_of_nodes); // Arguments to this function are the new z and x axes
  if (isnan(rot.r)) {
    rot = reb_rotation_identity();
  }
  struct reb_vec3d Omegarot = reb_vec3d_rotate(*Omega, rot);

  double spin_p;
  double theta_p;
  double phi_p;
  reb_tools_xyz_to_spherical(*Omega, &spin_p, &theta_p, &phi_p);
  double lr = laplace_radius(sim, rebx, planet, star);
  //printf("%f %f %f %f\n", planet->m / 3e-6, Omega->x, Omega->y, Omega->z);
  return theta_p - 0.5 * atan2(sin(2.*theta_p),cos(2.*theta_p) + (lr/d)*(lr/d)*(lr/d)*(lr/d)*(lr/d));
}

int main(int argc, char* argv[]){
    double ds[21] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};

    struct reb_simulation* sim = reb_simulation_create();
    sim->integrator         = REB_INTEGRATOR_WHFAST;
    sim->heartbeat          = heartbeat;

    ind = 0;
    if (argc == 2){
      strcat(title, argv[1]);
      strcat(title_remove, argv[1]);
      ind = atoi(argv[1]);
    }
    system(title_remove);

    // Load in the prescribed spin axis evolution
    int max_read = 100000; // so stupid
    double time[max_read];
    double spinx[max_read];
    double spiny[max_read];
    double spinz[max_read];

    FILE * ftime = fopen("140times.txt", "r");
    FILE * fsx = fopen("140sx.txt", "r");
    FILE * fsy = fopen("140sy.txt", "r");
    FILE * fsz = fopen("140sz.txt", "r");
    for (int i = 0; i < max_read; i++)
    {
        fscanf(ftime, "%lf", &time[i]);
        fscanf(fsx, "%lf", &spinx[i]);
        fscanf(fsy, "%lf", &spiny[i]);
        fscanf(fsz, "%lf", &spinz[i]);
    }

    // Create splines
    int n = 100000;
    sxfunc = rebx_create_interpolator(rebx, n, time, spinx, REBX_INTERPOLATION_SPLINE);
    syfunc = rebx_create_interpolator(rebx, n, time, spiny, REBX_INTERPOLATION_SPLINE);
    szfunc = rebx_create_interpolator(rebx, n, time, spinz, REBX_INTERPOLATION_SPLINE);

    // Initial conditions
    // Santerne et al 2019
    double mearth = 3e-6;
    double rho = 1.0 * pow(1.496e13, 3.) / (1.989e33); // 1 g/cm3 to rebound units
    struct reb_particle pf = {0};
    pf.m  = 12. * mearth;
    pf.r = pow(((3. * pf.m) / (4. * M_PI * rho)), 1./3.);
    reb_simulation_add(sim, pf);

    // Star
    double ms = 1.16;
    double af = 3.267176;
    reb_simulation_add_fmt(sim, "primary m a", pf, ms, af);

    rebx = rebx_attach(sim);
    //migration
    struct rebx_force* mof = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mof);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -1e7*2*M_PI*pow(af, -1.7));

    // tides spin
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    double planet_k2 = 0.6;
    struct reb_vec3d Omega_sv = {0};
    Omega_sv.x = spinx[0];
    Omega_sv.y = spiny[0];
    Omega_sv.z = spinz[0];
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.25 * pf.m * pf.r * pf.r);
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", Omega_sv);
    rebx_spin_initialize_ode(rebx, effect);

    // Laplace radius
    LR = laplace_radius(sim, rebx, &sim->particles[0], &sim->particles[1]);
    // Test particle
    double d = ds[ind] * LR;
    double le_theta = laplace_equilibrium(sim, rebx, &sim->particles[0], &sim->particles[1], d);

    //struct reb_orbit of = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    //struct reb_vec3d ofn = of.hvec;
    double mag_o;
    double theta_o;
    double phi_o;
    reb_tools_xyz_to_spherical(Omega_sv, &mag_o, &theta_o, &phi_o);

    reb_simulation_add_fmt(sim, "primary a inc Omega", pf, d, le_theta, phi_o + 90. * M_PI/180.);
    reb_simulation_move_to_com(sim);

    struct reb_orbit ot = reb_orbit_from_particle(sim->G, sim->particles[2], sim->particles[0]);
    struct reb_vec3d l = ot.hvec;
    double mag_l;
    double theta_l;
    double phi_l;
    reb_tools_xyz_to_spherical(l, &mag_l, &theta_l, &phi_l);

    struct reb_vec3d l_hat = reb_vec3d_normalize(l);
    struct reb_vec3d Omega_hat = reb_vec3d_normalize(Omega_sv);
    //printf("%f %f %f %f %f %f %f %f %f %f\n", l_hat.x, l_hat.y, l_hat.z, Omega_hat.x, Omega_hat.y, Omega_hat.z, theta_o * 180./M_PI, theta_l * 180./M_PI, phi_o * 180./M_PI, phi_l * 180./M_PI);
    //exit(1);

    tmax = 2e6 * 2 * M_PI;
    sim->dt = ot.P / 15.12345;
    sim->N_active=2;
    reb_simulation_integrate(sim, tmax);

    rebx_free_interpolator(sxfunc);
    rebx_free_interpolator(syfunc);
    rebx_free_interpolator(szfunc);
    rebx_free(rebx);

    reb_simulation_free(sim);
}

void heartbeat(struct reb_simulation* sim){

    struct reb_vec3d new_Omega = {0};
    new_Omega.x = rebx_interpolate(rebx, sxfunc, sim->t);
    new_Omega.y = rebx_interpolate(rebx, syfunc, sim->t);
    new_Omega.z = rebx_interpolate(rebx, szfunc, sim->t);
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", new_Omega);

    // Output spin and orbital information to file
    if(reb_simulation_output_check(sim, 100. * 2 * M_PI)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      struct reb_particle* planet = &sim->particles[0];
      struct reb_particle* sun = &sim->particles[1];
      struct reb_particle* test = &sim->particles[2];

      struct reb_orbit of = reb_orbit_from_particle(sim->G, *sun, *planet);
      struct reb_vec3d nf = of.hvec;

      struct reb_orbit ot = reb_orbit_from_particle(sim->G, *test, *planet);
      struct reb_vec3d nt = ot.hvec;
      double at = ot.a;

      struct reb_vec3d* Omega_p_inv = rebx_get_param(rebx, planet->ap, "Omega");
      struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, nf);
      struct reb_rotation rot = reb_rotation_init_to_new_axes(nf, line_of_nodes); // Arguments to this function are the new z and x axes
      if (isnan(rot.r)) {
        rot = reb_rotation_identity();
      }
      struct reb_vec3d nt_rot = reb_vec3d_rotate(nt, rot);

      double mag_p;
      double theta_p;
      double phi_p;
      reb_tools_xyz_to_spherical(nt_rot, &mag_p, &theta_p, &phi_p);

      FILE* sf = fopen(title, "a");
      fprintf(sf, "%f,%f,%f,%f,%e,%e,%e\n",sim->t,Omega_p_inv->x, Omega_p_inv->y, Omega_p_inv->z,theta_p,phi_p,at/planet->r);
      fclose(sf);


    }


    //if(reb_simulation_output_check(sim, 10.)){        // outputs to the screen
    //    reb_simulation_output_timing(sim, tmax);
    //}

    struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    if (o.a < 1.37 && first_check){
      rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
      first_check=0;
    }
}
