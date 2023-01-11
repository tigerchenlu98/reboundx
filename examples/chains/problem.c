/**
 * Self-Consistent Spin-Tidal and Dynamical equations of motion (Eggleton et. al 1998)
 *
 * This example shows how to add quadrupole and tidal distortions to bodies with structure, letting us consistently track the spin-axis and dynamical evolution of the system.
 * In particular, this simulates the pseudo-synchronization of a fiducial hot Jupiter.
 * The hot Jupiter is initialized with a slight eccentricity, nontrivial obliquity and fast rotation.
 * Under the influence of tidal dissipation, we see the following rapidly occur: circularization of the orbit, obliquity damping down to 0, and rotation settling to the pseudo-synchronous value predicted by Hut (1981)
 * Definitely see the corresponding ipython example, as well as the documentation, for more in-depth explanations regarding the various parameters that may be set in this simulation.
 * Also see Lu et. al (in review), Eggleton et. al (1998), Hut (1981).
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "tides_spin.c"

void heartbeat(struct reb_simulation* sim);
double tmax = 1e5 * 2 * M_PI;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();

    // Star: M-dwarf
    const double solar_mass = 0.08;
    const double solar_rad = 0.117 * 0.00465; // Bodies with structure require radius! This is one solar radius
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    // Trappist-1g
    const double g_mass = 1.321 * 3e-6;
    const double g_rad = 1.129 * 4.264e-5;
    const double g_a = 0.04508048 * 1.09;
    const double g_e = 0.061;
    const double g_inc = 0.37 * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", g_mass, g_rad, g_a, g_e, g_inc); // Planet 1

    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    double ts = orb.P / 15.; // timestep as a function of orbital period of innermost planet

    // Trappist-1h
    const double h_mass = 0.326 * 3e-6;
    const double h_rad = 0.755 * 4.264e-5;
    const double h_a = 0.06189 * 1.1;
    const double h_e = 0.0;
    const double h_inc = 0.195 * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", h_mass, h_rad, h_a, h_e, h_inc);

    sim->N_active = 3;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = ts;
    sim->heartbeat = heartbeat;

    // Add tides_spin as a REBOUNDx additional effect
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    // Star
    const double solar_k2 = 0.307;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);

    const double solar_spin_period = 3.3 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin});

    const double rg_star = 0.2;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", rg_star * solar_mass * solar_rad * solar_rad);
    //rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

    // Planets - all assumed to have same tidal parameters
    const double planet_k2 = 0.299;
    const double planet_tau = (712.37 / 86400) * 2 * M_PI / 365; // seconds->days->reb years
    const double planet_rg = 0.3308;

    for (int i = 1; i < sim->N_active; i++){
      rebx_set_param_double(rebx, &sim->particles[i].ap, "k2", planet_k2);
      rebx_set_param_double(rebx, &sim->particles[i].ap, "tau", planet_tau);
      rebx_set_param_double(rebx, &sim->particles[i].ap, "I", planet_rg * sim->particles[i].m * sim->particles[i].r * sim->particles[i].r);

      // tidally locked & 0 obliquity
      struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[i], sim->particles[0]);
      const double spin_period = orb.P;
      const double spin_rate = 2 * M_PI / spin_period;

      struct reb_vec3d norm = orb.hvec;
      struct reb_vec3d nhat = reb_vec3d_normalize(norm);
      struct reb_vec3d Omega_p = reb_vec3d_mul(nhat, spin_rate);
      rebx_set_param_vec3d(rebx, &sim->particles[i].ap, "Omega", Omega_p);
    }

    // Migration
    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -3e5 * 2 * M_PI);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", (-3e5 * 2 * M_PI) / 1.1);

    reb_move_to_com(sim);

    struct reb_vec3d newz = rebx_tools_total_angular_momentum(rebx);
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);

    rebx_simulation_irotate(rebx, rot);
    rebx_spin_initialize_ode(rebx, effect);

    system("rm -v output.txt"); // remove previous output file
    FILE* of = fopen("output.txt", "a");
    fprintf(of, "t,a1,i1,e1,pom1,Om1,mag1,theta1,phi1,a2,i2,e2,pom2,Om2,mag2,theta2,phi2\n");
    fclose(of);
    reb_integrate(sim, tmax/2);

    printf("\nMigration Switching Off\n");
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", INFINITY);

    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 10)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("output.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      struct reb_particle* star = &sim->particles[0];
      fprintf(of, "%e", sim->t);
      for (int i = 1; i < sim->N_active; i++){
        struct reb_particle* p = &sim->particles[i];

        struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, *p, *star);
        double a = orb.a;
        double Om = orb.Omega;
        double inc = orb.inc;
        double pom = orb.pomega;
        double e = orb.e;

        // The spin vector in the inertial (in this case, invariant frame)
        struct reb_vec3d* Omega_inv = rebx_get_param(rebx, p->ap, "Omega");

        // Transform spin vector into planet frame, w/ z-axis aligned with orbit normal and x-axis aligned with line of nodes
        struct reb_vec3d orbit_normal = orb.hvec;
        struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, orbit_normal);
        struct reb_rotation rot = reb_rotation_init_to_new_axes(orbit_normal, line_of_nodes); // Arguments to this function are the new z and x axes
        struct reb_vec3d srot = reb_vec3d_rotate(*Omega_inv, rot); // spin vector in the planet's frame


        // Interpret the spin axis in the more natural spherical coordinates
        double mag;
        double theta;
        double phi;
        reb_tools_xyz_to_spherical(srot, &mag, &theta, &phi);
        fprintf(of, ",%e,%e,%e,%e,%e,%e,%e,%e", a, inc, e, pom, Om, mag, theta, phi);
      }
      fprintf(of, "\n");
      fclose(of);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
