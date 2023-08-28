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
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Initial conditions
    // Setup constants
    sim->integrator         = REB_INTEGRATOR_IAS15; // IAS15 is used for its adaptive timestep:
                                                   // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                   // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
    sim->heartbeat          = heartbeat;
    sim->N_active = 3;
    //sim->collision = REB_COLLISION_NONE;

    // Initial conditions
    // Planet
    struct reb_particle star = {0};
    star.m  = 1.0;
    star.r = 1.0 * 0.00465;
    reb_add(sim, star);

    // Planet
    double pm = 5. * 3e-6;
    double pa1 = 0.55;
    double pa2 = pa1 * pow(1.55, 2./3.);
    double pe = 0.01;
    double pr = 1.5*2.5*4.26e-5;
    double pinc1 = 0.5 * (M_PI / 180.);
    double pinc2 = -0.5 * (M_PI / 180.);
    reb_add_fmt(sim, "m a r e inc", pm, pa1, pr, pe, pinc1);
    reb_add_fmt(sim, "m a r e inc", pm, pa2, pr, pe, pinc2);



    struct reb_particle planet = sim->particles[1];
    //double testa = 1.6 * pr;
    //reb_add_fmt(sim, "primary a", planet, testa);

    struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, sim->particles[2], sim->particles[0]);

    // timestep
    //sim->dt = o1.P / 15.12345;


    // Add REBOUNDx effects
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    struct rebx_force* mof = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mof);

    double tau1 = -2.5e7 * (2 * M_PI);
    double tau2 = tau1 / 1.1;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", tau1);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", tau2);

    // Sun
    const double solar_k2 = 0.07;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * star.m * star.r * star.r);
    //rebx_set_param_double(rebx, &sim->particles[0].ap, "tau_a", -1e5);

    const double solar_spin_period = 20. * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default

    // Planet
    const double p_c = 0.5;
    const double p_k2 = 0.5;
    const double pQ = 5e2;
    const double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = 2 * M_PI / spin_period_p;
    const double theta_p = 1.0 * M_PI / 180.;
    const double phi_p = 0. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);

    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", p_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", p_c * pm * pr * pr);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1/(2. * pQ * o1.n));

    rebx_set_param_vec3d(rebx, &sim->particles[2].ap, "Omega", Omega_sv);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", p_k2);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "I", p_c * pm * pr * pr);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau", 1/(2. * pQ * o2.n));

    // Damping
    //rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_i", 1e3 * 2 * M_PI);


    reb_move_to_com(sim);

    // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
    // nodes with the invariable plane (along z cross newz)
    struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    if (isnan(rot.r)) {
      rot = reb_rotation_identity();
    }
    rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)

    rebx_spin_initialize_ode(rebx, effect);


    system("rm -v test.txt");        // delete previous output file
    //system("rm -v evol.txt");
    FILE* of = fopen("test.txt", "w");
    fprintf(of, "t,a1,a2,Omega1,Omega2,mag_p1,theta_p1,mag_p2,theta_p2\n");
    //for (int i = sim->N_active; i < sim->N; i++){
    //  fprintf(of, ",tx%d,ty%d,tz%d,a%d,e%d,hx%d,hy%d,hz%d",i-1,i-1,i-1,i-1,i-1,i-1,i-1,i-1);
    //}
    //fprintf(of, "\n");
    fclose(of);

    //struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    tmax = 4e6*2*M_PI;
    reb_integrate(sim, 3./2.*tmax);

    // turn off migration
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", INFINITY);
    reb_integrate(sim, 10.*tmax);

    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 10. * M_PI)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("test.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* p1 = &sim->particles[1];
      struct reb_particle* p2 = &sim->particles[2];
      struct reb_particle* t = &sim->particles[3];

      struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
      struct reb_vec3d* Omega_p1 = rebx_get_param(rebx, p1->ap, "Omega");
      struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *p2, *sun);
      struct reb_vec3d* Omega_p2 = rebx_get_param(rebx, p2->ap, "Omega");
      struct reb_orbit ot = reb_tools_particle_to_orbit(sim->G, *t, *p1);

      struct reb_vec3d line_of_nodes_1 = reb_vec3d_cross((struct reb_vec3d){.z =1}, o1.hvec);
      struct reb_rotation rot1 = reb_rotation_init_to_new_axes(o1.hvec, line_of_nodes_1); // Arguments to this function are the new z and x axes
      struct reb_vec3d srot1 = reb_vec3d_rotate(*Omega_p1, rot1); // spin vector in the planet's frame

      // Interpret the spin axis in the more natural spherical coordinates
      double mag_p1;
      double theta_p1;
      double phi_p1;
      reb_tools_xyz_to_spherical(srot1, &mag_p1, &theta_p1, &phi_p1);

      struct reb_vec3d line_of_nodes_2 = reb_vec3d_cross((struct reb_vec3d){.z =1}, o2.hvec);
      struct reb_rotation rot2 = reb_rotation_init_to_new_axes(o2.hvec, line_of_nodes_2); // Arguments to this function are the new z and x axes
      struct reb_vec3d srot2 = reb_vec3d_rotate(*Omega_p2, rot2); // spin vector in the planet's frame

      // Interpret the spin axis in the more natural spherical coordinates
      double mag_p2;
      double theta_p2;
      double phi_p2;
      reb_tools_xyz_to_spherical(srot2, &mag_p2, &theta_p2, &phi_p2);

      fprintf(of, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",sim->t,o1.a,o2.a,o1.Omega,o2.Omega,mag_p1,theta_p1,mag_p2,theta_p2,ot.x,ot.y,ot.z);
      //for (int i = sim->N_active; i < sim->N; i++){
      //    struct reb_particle* test = &sim->particles[i];
      //    struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *test, *p);
      //    fprintf(of, ",%f,%f,%f,%f,%f,%f,%f,%f",test->x,test->y,test->z,o.a,o.e,o.hvec.x, o.hvec.y,o.hvec.z);
      //}
      //fprintf(of, "\n");
      fclose(of);
      /*
      FILE* ef = fopen("evol.txt", "a");
      struct reb_particle* test = &sim->particles[2];
      double dx = p->x - test->x;
      double dy = p->y - test->y;
      double dz = p->z - test->z;
      double d = sqrt(dx*dx+dy*dy+dz*dz);
      double vx = p->vx - test->vx;
      double vy = p->vy - test->vy;
      double vz = p->vz - test->vz;
      double v = sqrt(vx*vx+vy*vy+vz*vz);
      double ax = p->ax - test->ax;
      double ay = p->ay - test->ay;
      double az = p->az - test->az;
      double a = sqrt(ax*ax+ay*ay+az*az);
      fprintf(ef, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", sim->t,dx,dy,dz,d,vx,vy,vz,v,ax,ay,az,a);
      fclose(ef);
      */
    }

    //if(reb_output_check(sim, 0.1*M_PI)){        // outputs to the screen
    //    reb_output_timing(sim, tmax);
    //}
}
