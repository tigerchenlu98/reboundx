/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 *
 * This example includes self-consistent spin, tidal & dynamical effects
 * as well as general relativity
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax = 1e5 * 2 * M_PI;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants
    sim->dt             = 2*M_PI*1.;     // initial timestep
    sim->integrator        = REB_INTEGRATOR_IAS15;
    //sim->heartbeat        = heartbeat;

    // Initial conditions

    struct reb_particle star = {0};
    star.m  = 1.05;
    star.r = 0.968 * 0.00465;
    reb_add(sim, star);

    double planet_m  = 4.2 * 9.55e-4; // in Jupiter masses
    double planet_r = 0.974 * 4.676e-4;
    double planet_a = 5.0;
    double planet_e = 0.1;
    reb_add_fmt(sim, "m r a e", planet_m, planet_r, planet_a, planet_e);

    // The perturber
    double perturber_inc = 85.6 * (M_PI / 180.);
    double perturber_mass = 1.1;
    double perturber_a  = 1000.;
    double perturber_e = 0.5;
    reb_add_fmt(sim, "m a e inc", perturber_mass, perturber_a, perturber_e, perturber_inc);

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);
    // Sun
    const double solar_spin_period = 20 * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_k2 = 0.028;
    const double solar_tau = (0.2 / solar_k2) * (1 / 86400.) * 2 * M_PI / 365.; // seconds to reb years
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.08 * star.m * star.r * star.r);
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin});

    //double solar_Q = 1e6;
    //struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    //double solar_tau = 1 / (2 * solar_Q * orb.n);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);
    printf("solar tau: %e\n", solar_tau);

    // P1
    const double spin_period_p = (10./24.) * 2. * M_PI / 365.; // days to reb years
    const double spin_p = (2. * M_PI) / spin_period_p;
    const double planet_k2 = 0.51;
    const double planet_Q = 3e6;
    const double planet_tau = (0.02 / planet_k2) * (1 / 86400.) * 2 * M_PI / 365.; // seconds to reb years
    const double theta_1 = 0. * M_PI / 180.;
    const double phi_1 = 0. * M_PI / 180;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * planet_m * planet_r * planet_r);
    
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_1, phi_1);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);
    //double planet_tau = 1. / (2. * planet_Q * orb.n);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", planet_tau);
   
    printf("planet tau: %e\n", planet_tau);


    // add GR:
    struct rebx_force* gr = rebx_load_force(rebx, "gr_full");
    rebx_add_force(rebx, gr);

    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units

    reb_move_to_com(sim);
    
    struct reb_vec3d newz = rebx_tools_total_angular_momentum(rebx);
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z=1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz,newx);
    rebx_simulation_irotate(rebx, rot);

    rebx_spin_initialize_ode(rebx, effect);

    FILE* f = fopen("01_10_winn09_lec10_biggest_tau.txt","w");
    fprintf(f, "t,star_sx,star_sy,star_sz,magstar,a1,i1,e1,s1x,s1y,s1z,mag1,pom1,Om1,f1,a2,i2,e2,Om2,pom2,o1x,o1y,o1z,o2x,o2y,o2z\n");

    for (int i=0; i<10000000; i++){
        struct reb_particle* sun = &sim->particles[0];
        struct reb_particle* p1 = &sim->particles[1];
        struct reb_particle* pert = &sim->particles[2];

        struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");
        double star_sx = Omega_sun->x;
        double star_sy = Omega_sun->y;
        double star_sz = Omega_sun->z;
        double magstar = sqrt(star_sx * star_sx + star_sy * star_sy + star_sz * star_sz);

        struct reb_vec3d* Omega_p = rebx_get_param(rebx, p1->ap, "Omega");
        double s1x = Omega_p->x;
        double s1y = Omega_p->y;
        double s1z = Omega_p->z;
        double mag1 = sqrt(s1x * s1x +s1y * s1y + s1z * s1z);
            

        struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
        double a1 = o1.a;
        double Om1 = o1.Omega;
        double i1 = o1.inc;
        double pom1 = o1.pomega;
        double e1 = o1.e;
        double f1 = o1.f;
        struct reb_vec3d norm1 = o1.hvec;

        struct reb_particle com = reb_get_com_of_pair(sim->particles[0],sim->particles[1]);
        struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *pert, com);
        double a2 = o2.a;
        double Om2 = o2.Omega;
        double i2 = o2.inc;
        double pom2 = o2.pomega;
        double e2 = o2.e;
        struct reb_vec3d norm2 = o2.hvec;

        fprintf(f, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%.e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, star_sx, star_sy, star_sz, magstar, a1, i1, e1, s1x, s1y, s1z, mag1, pom1, Om1, f1, a2, i2, e2, Om2, pom2, norm1.x, norm1.y, norm1.z, norm2.x, norm2.y, norm2.z);
        reb_integrate(sim, sim->t+(100 * 2 * M_PI));
    }
   rebx_free(rebx);
   reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        // reb_output_timing(r, tmax);
    }
    /*
    if(reb_output_check(r, 12.)){            // outputs to a file
        reb_output_orbits(r, "orbits.txt");
    }
    */
}
