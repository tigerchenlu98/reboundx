/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 *
 * This is the same Kozai example implemented in Lu et. al (2023)
 * Also, see the ipython examples prefixed TidesSpin for in-depth exploration of the parameters that can be set in this simulation.
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax = 1e4 * M_PI * 2.;
int Ntest = 10;
double alignment_ts = (1e3 * M_PI * 2.);
double inv_alignment_ts = 1./(1e100 * M_PI * 2.);

void derivatives(struct reb_ode* const ode, double* const yDot, const double* const y, const double t){
    // From Zanazzi & Lai 2017
    // hard-coded relevant physical parameters
    const double ms = ode->r->particles[0].m;
    const double mp = ode->r->particles[1].m;
    const double rp = ode->r->particles[1].r;
    struct reb_orbit op = reb_tools_particle_to_orbit(ode->r->G, ode->r->particles[1], ode->r->particles[0]); // planet orbit
    const double opa = op.a;
    const double ope = op.e;
    struct reb_vec3d lp = op.hvec;
    struct reb_vec3d lp_hat = reb_vec3d_normalize(lp);

    const struct reb_vec3d* const sp_ptr = rebx_get_param(ode->r->extras, ode->r->particles[1].ap, "Omega");
    struct reb_vec3d sp;
    if (sp_ptr != NULL){
      sp = *sp_ptr;
    }
    struct reb_vec3d sp_hat = reb_vec3d_normalize(sp);

    // This stuff for alignment torque
    struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, op.hvec);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(op.hvec, line_of_nodes); // Invariant plane to planet
    if (isnan(rot.r)) {
      rot = reb_rotation_identity();
    }
    struct reb_vec3d sp_rot = reb_vec3d_rotate(sp, rot);

    double omega;
    double theta_p;
    double phi_p;
    reb_tools_xyz_to_spherical(sp_rot, &omega, &theta_p, &phi_p); // theta and phi in the planet frame
    struct reb_rotation invrot = reb_rotation_inverse(rot); // planet frame to invariant plane

    // Planet quadrupole
    const double* k2 = rebx_get_param(ode->r->extras, ode->r->particles[1].ap, "k2");
    const double j2 = (*k2 / 3.) * (omega * omega * rp * rp * rp) / (ode->r->G * mp);
    const double lr = pow(2. * j2 * rp * rp * opa * opa * opa * pow((1 - ope * ope),(3./2.)) * mp / ms, (1./5.));

    for (unsigned int i = 0; i < Ntest; i++){
      // Multiple particles
      struct reb_vec3d l_hat = {};
      l_hat.x = y[i*4+0];
      l_hat.y = y[i*4+1];
      l_hat.z = y[i*4+2];

      // Planet-TP distance
      double d = y[i*4+3];

      // star
      const double star_prefactor = ((3. * ode->r->G * ms * d * d) / (4 * opa * opa * opa)) * reb_vec3d_dot(l_hat, lp_hat);;
      struct reb_vec3d l_cross_lp = reb_vec3d_cross(l_hat, lp_hat);
      struct reb_vec3d tstar = reb_vec3d_mul(l_cross_lp, star_prefactor);

      // planet
      const double planet_prefactor = ((3. * ode->r->G * mp * rp * rp * j2) / (2. * d * d * d)) * reb_vec3d_dot(l_hat, sp_hat);
      struct reb_vec3d l_cross_sp = reb_vec3d_cross(l_hat, sp_hat);
      struct reb_vec3d tplanet = reb_vec3d_mul(l_cross_sp, planet_prefactor);

      // simple alignment torque
      // Align towards Laplace Equilibrium in planet frame

      const double le_theta = theta_p - 0.5 * atan2(sin(2.*theta_p),cos(2.*theta_p) + (lr/d)*(lr/d)*(lr/d)*(lr/d)*(lr/d)); // In the planet frame
      struct reb_vec3d beta_vec = reb_tools_spherical_to_xyz(1., le_theta, phi_p);// In the planet frame
      reb_vec3d_irotate(&beta_vec, invrot); // rotates into xyz frame.

      struct reb_vec3d delta_i = {};
      delta_i.x = beta_vec.x - l_hat.x;
      delta_i.y = beta_vec.y - l_hat.y;
      delta_i.z = beta_vec.z - l_hat.z;

      const double align_prefactor = inv_alignment_ts;// * ode->r->dt;//inv_alignment_ts * exp(-1. * ode->r->t * inv_alignment_ts);

      struct reb_vec3d talign = {};//reb_vec3d_mul(delta_i, align_prefactor);//inv_alignment_ts);

      // DiffEq
      const double inv_prefactor = 1. / (d * d * sqrt(ode->r->G * mp / (d * d * d)));
      yDot[i*4+0] = inv_prefactor * (tstar.x + tplanet.x) + talign.x;
      yDot[i*4+1] = inv_prefactor * (tstar.y + tplanet.y) + talign.y;
      yDot[i*4+2] = inv_prefactor * (tstar.z + tplanet.z) + talign.z;
      yDot[i*4+3] = 0;

    }
}


//char title[100] = "831_fd_";
//double semis[] = {1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0};

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Initial conditions
    // Setup constants
    sim->integrator         = REB_INTEGRATOR_WHFAST; // IAS15 is used for its adaptive timestep:
                                                   // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                   // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
    sim->heartbeat          = heartbeat;
    sim->N_active = 2;
/*
    int index = 0;
    if (argc == 2){
       strcat(title, argv[1]);
       index = atoi(argv[1]);
    }
    double val = semis[index];
*/
    // Initial conditions
    // Sun
    struct reb_particle star = {0};
    double star_m = 1.0;
    double star_r = 1.0 * 0.00465;
    star.m  = star_m;
    star.r = star_r;
    reb_add(sim, star);

    // Planet
    //struct reb_particle planet = {0};
    double pm = 5 * 3e-6;
    double pa = 0.5;
    double pr = 2.5 * 4.26e-5;
    //planet.m = pm;
    //planet.r = pr;
    reb_add_fmt(sim, "m a r", pm, pa, pr);

    struct reb_particle planet = sim->particles[1];

// ----------------HARMONIC OSCILLATOR----------------------
    struct reb_ode* ho = reb_create_ode(sim,Ntest*4+1);   // Add an ODE with 2 dimensions
    ho->derivatives = derivatives;              // Right hand side of the ODE
// ----------------HARMONIC OSCILLATOR----------------------

    double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = 2 * M_PI / spin_period_p;
    const double theta_p = 45. * M_PI / 180.;

    const double k2 = 0.4;

    const double j2 = (k2 / 3.) * (spin_p * spin_p * pr * pr * pr) / (sim->G * pm);
    const double lr = pow(2. * j2 * pr * pr * pa * pa * pa * pow((1 - 0.0),(3./2.)) * pm / star_m, (1./5.));

    // Test particle
    for (unsigned int i = 0; i < Ntest; i++){
      double d = (1.5 + i) * pr;
      double le_theta = theta_p - 0.5 * atan2(sin(2.*theta_p),cos(2.*theta_p) + (lr/d)*(lr/d)*(lr/d)*(lr/d)*(lr/d));
      reb_add_fmt(sim, "primary a inc Omega", planet, d, le_theta, 90. * M_PI/180.);

      struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[2], sim->particles[1]);
      struct reb_vec3d l = o.hvec;
      struct reb_vec3d l_hat = reb_vec3d_normalize(l);

      //if (i == 0){
      //  sim->dt = o.P / 10.6789;
      //}

      ho->y[i*4+0] = l_hat.x;                               // Initial conditions
      ho->y[i*4+1] = l_hat.y;
      ho->y[i*4+2] = l_hat.z;
      ho->y[i*4+3] = d;

      //printf("%f %f %f\n", l_hat.x, l_hat.y, l_hat.z);

      // remove particle
      reb_remove(sim, 2, 1);
    }

    struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    sim->dt = o.P / 15.6789;



    // Add REBOUNDx effects
    // First tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    // Planet
    const double phi_p = 0. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);

    //struct reb_vec3d Omega_align = reb_tools_spherical_to_xyz(spin_p, 45. * M_PI/180., phi_p);
    //rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "alignment", Omega_align);
    //rebx_set_param_double(rebx, &sim->particles[1].ap, "alignment_ts", 1e4*2*M_PI);

    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.5 * pm * pr * pr);

    reb_move_to_com(sim);
    rebx_spin_initialize_ode(rebx, effect);

    system("rm -v output_ho_align.txt");        // delete previous output file
    FILE* of = fopen("output_ho_align.txt", "w");
    fprintf(of, "t,px,py,pz,omega_p,theta_p,phi_p");
    //fprintf(of, "t,px,py,pz");
    for (int i = 1; i < Ntest+1; i++){
      //fprintf(of, ",nx%d,ny%d,nz%d,a%d,e%d,theta%d,phi%d",i,i,i,i,i,i,i);
      fprintf(of, ",nx%d,ny%d,nz%d,a%d,theta%d,phi%d",i,i,i,i,i,i);
    }
    fprintf(of,",lr\n");
    fclose(of);

    reb_integrate(sim, tmax);

    //reb_free_ode(ho);
    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 1. * 2 * M_PI)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("output_ho_align.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      const struct reb_vec3d* const sp_ptr = rebx_get_param(sim->extras, sim->particles[1].ap, "Omega");
      struct reb_vec3d sp;
      if (sp_ptr != NULL){
        sp = *sp_ptr;
      }
      struct reb_vec3d sp_hat = reb_vec3d_normalize(sp);


      struct reb_ode** odes = sim->odes;
      struct reb_ode* lhat_ode = odes[0];



      struct reb_orbit op = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]); // planet orbit
      // for planet frame
      //struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, op.hvec);

      // for aligned with spin axis
      /*
      struct reb_vec3d newx = {};
      newx.x = sp_hat.x;
      newx.y = sp_hat.y;
      */
/*
      struct reb_rotation rot = reb_rotation_init_to_new_axes(op.hvec, newx); // Invariant plane to planet
      if (isnan(rot.r)) {
        rot = reb_rotation_identity();
      }
      reb_vec3d_irotate(&sp_hat, rot);
      reb_vec3d_irotate(&sp, rot);
      */

      double omega;
      double theta_p;
      double phi_p;
      reb_tools_xyz_to_spherical(sp, &omega, &theta_p, &phi_p);

      fprintf(of, "%f,%f,%f,%f,%f,%f,%f",sim->t,sp_hat.x, sp_hat.y, sp_hat.z,omega,theta_p,phi_p);
      //fprintf(of, "%f,%f,%f,%f",sim->t,sp_hat.x, sp_hat.y, sp_hat.z);

      for (unsigned int i = 0; i < Ntest; i++){

        struct reb_vec3d lhat = {};
        lhat.x = lhat_ode->y[i*4+0];
        lhat.y = lhat_ode->y[i*4+1];
        lhat.z = lhat_ode->y[i*4+2];

        double mag;
        double theta;
        double phi;

        //reb_vec3d_irotate(&lhat, rot); // rotates into planet frame
        reb_tools_xyz_to_spherical(lhat, &mag, &theta, &phi);

        fprintf(of, ",%e,%e,%e,%e,%f,%f", lhat.x, lhat.y, lhat.z, lhat_ode->y[i*4+3],theta,phi);
        //fprintf(of, ",%e,%e,%e,%e,%f,%f,%f", lhat.x, lhat.y, lhat.z, ot.a, ot.e, theta,phi);
      }


      const double* k2 = rebx_get_param(sim->extras, sim->particles[1].ap, "k2");
      double rp = sim->particles[1].r;
      double mp = sim->particles[1].m;
      double ms = sim->particles[0].m;
      const double j2 = ((*k2) / 3.) * (omega * omega * rp * rp * rp) / (sim->G * mp);
      const double lr = pow(2. * j2 * rp * rp * op.a * op.a * op.a * pow((1 - op.e * op.e),(3./2.)) * mp / ms, (1./5.));
      //const double lr = pow(2. * j2 * (mp / ms) * rp * rp * op.a * op.a * op.a, 1./5.);

      fprintf(of, ",%f\n", lr);
      fclose(of);
      //zexit(1);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
