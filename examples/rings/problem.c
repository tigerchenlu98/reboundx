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
    sim->integrator         = REB_INTEGRATOR_WHFAST; // IAS15 is used for its adaptive timestep:
                                                   // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                   // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
    sim->heartbeat          = heartbeat;
    sim->N_active = 2;
    sim->collision = REB_COLLISION_NONE;

    // Initial conditions
    // Planet
    struct reb_particle planet = {0};
    planet.m  = 10. * 3e-6;
    planet.r = 3. * 4.26e-5;
    reb_add(sim, planet);

    // Sun
    double pm = 1.0;
    double pa = 0.5;
    double pe = 0.0;
    double pr = 1.0 * 0.00465;
    double pinc = 0. * (M_PI / 180.);
    reb_add_fmt(sim, "m a r e inc", pm, pa, pr, pe, pinc);


    // Test particle
    double inner_limit = 1.5 * planet.r;
    double rin = 2.0 * planet.r;
    double rout = 5.0 * planet.r;//0.2 * pa * pow((pm / (3 * star.m)), (1./3.)); //hill radius

    const double theta_p = 0.0 * M_PI / 180.;
    int Ntest = 12;


    for (int i = 0; i < Ntest; i++){
      double testa = rin + i * ((rout - rin) / Ntest);
      reb_add_fmt(sim, "primary a", planet, testa);
    }

    // timestep
    double period = sqrt(((4 * M_PI * M_PI) / planet.m) * inner_limit * inner_limit * inner_limit);
    sim->dt = period / 10.12345;


    // Add REBOUNDx effects
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    struct rebx_force* damping = rebx_load_force(rebx, "laplace_damping");
    rebx_add_force(rebx, effect);
    //rebx_add_force(rebx, damping);

    //struct rebx_force* mof = rebx_load_force(rebx, "modify_orbits_forces");
    //rebx_add_force(rebx, mof);

    // Sun
    const double solar_k2 = 0.01;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.07 * pm * pr * pr);

    const double solar_spin_period = 27. * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default

    // Planet
    double spin_period_p = 5. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = 2 * M_PI / spin_period_p;
    const double phi_p = 0. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", Omega_sv);

    const double p_k2 = 0.5;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", p_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.225 * planet.m * planet.r * planet.r);

    // Laplace Radius
    double j2 = (spin_p * spin_p * planet.r * planet.r * planet.r) / (3. * planet.m) * (p_k2);
    double laplace_radius = pow(j2 * planet.r * planet.r * pa*pa*pa * pow((1 - pe*pe),(3./2.)) * planet.m / pm, (1./5.));
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lr", laplace_radius);
    rebx_set_param_int(rebx, &sim->particles[0].ap, "primary", 1);

    // Planetary Spin alignment torque
    struct reb_vec3d alignment_vec = reb_tools_spherical_to_xyz(1., 45. * M_PI / 180., 0.);
    const double alignment_ts = 1e5;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "alignment", alignment_vec);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "alignment_ts", alignment_ts);


    // Damping
    for (int i = 0; i < Ntest; i++){
      rebx_set_param_double(rebx, &sim->particles[i+2].ap, "tau_i", 1e3 * 2. * M_PI);
    }

    rebx_set_param_double(rebx, &damping->ap, "ide_position", 2.2*planet.r);
    rebx_set_param_double(rebx, &damping->ap, "ide_width", 0.2*planet.r);

    //rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_inc", 10. * 2 * M_PI);


    reb_move_to_com(sim);
    rebx_spin_initialize_ode(rebx, effect);


    //system("rm -v output.txt");        // delete previous output file
    //system("rm -v evol.txt");
    FILE* of = fopen("output_fast_damp.txt", "w");
    fprintf(of, "t,px,py,pz,Omx,Omy,Omz,nx,ny,nz");
    for (int i = sim->N_active; i < sim->N; i++){
      fprintf(of, ",tx%d,ty%d,tz%d,a%d,e%d,hx%d,hy%d,hz%d",i-1,i-1,i-1,i-1,i-1,i-1,i-1,i-1);
    }
    fprintf(of, "\n");
    fclose(of);

    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    tmax = 1e6;
    reb_integrate(sim, tmax);

    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 100.)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("output_fast_damp.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      struct reb_particle* p = &sim->particles[0];
      struct reb_particle* sun = &sim->particles[1];
      struct reb_vec3d* Omega_p = rebx_get_param(rebx, p->ap, "Omega");
      struct reb_orbit op = reb_tools_particle_to_orbit(sim->G, *p, *sun);
      fprintf(of, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",sim->t,p->x,p->y,p->z,Omega_p->x,Omega_p->y,Omega_p->z,op.hvec.x,op.hvec.y,op.hvec.z);
        for (int i = sim->N_active; i < sim->N; i++){
          struct reb_particle* test = &sim->particles[i];
          struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *test, *p);
          fprintf(of, ",%f,%f,%f,%f,%f,%f,%f,%f",test->x,test->y,test->z,o.a,o.e,o.hvec.x, o.hvec.y,o.hvec.z);
      }
      fprintf(of, "\n");
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

    if(reb_output_check(sim, 0.1*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
