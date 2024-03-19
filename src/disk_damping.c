#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static struct reb_vec3d rebx_calculate_disk_damping(struct reb_simulation* const sim, struct reb_particle* p, struct reb_particle* planet, struct reb_particle* primary, double tau_i, double k2, double mag_p, double theta_p, double phi_p, struct reb_rotation invrot){

    //struct reb_orbit op = reb_tools_particle_to_orbit_err(sim->G, *planet, *sun, &err);
    // Planet orbit
    // semimajor axis
    const double mu_p = sim->G*(planet->m+primary->m);
    const double dx_p = planet->x - primary->x;
    const double dy_p = planet->y - primary->y;
    const double dz_p = planet->z - primary->z;
    const double dvx_p = planet->vx - primary->vx;
    const double dvy_p = planet->vy - primary->vy;
    const double dvz_p = planet->vz - primary->vz;
    const double d_p = sqrt ( dx_p*dx_p + dy_p*dy_p + dz_p*dz_p );
    const double v2_p = dvx_p*dvx_p + dvy_p*dvy_p + dvz_p*dvz_p;
    const double vcirc2_p = mu_p/d_p;
    const double diff2_p = v2_p - vcirc2_p;
    const double ap = -mu_p/( v2_p - 2.*vcirc2_p );	// semi major axis

    // eccentricity
    const double vr_p = (dx_p*dvx_p + dy_p*dvy_p + dz_p*dvz_p)/d_p;
    const double rvr_p = d_p*vr_p;
    const double muinv_p = 1./mu_p;

    const double ex = muinv_p*( diff2_p*dx_p - rvr_p*dvx_p );
    const double ey = muinv_p*( diff2_p*dy_p - rvr_p*dvy_p );
    const double ez = muinv_p*( diff2_p*dz_p - rvr_p*dvz_p );
    const double ep = sqrt(ex*ex+ey*ey+ez*ez);

    // Test particle orbit
    //struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *source, &err);
    const double mu = sim->G*(planet->m+p->m);
    const double dx = p->x - planet->x;
    const double dy = p->y - planet->y;
    const double dz = p->z - planet->z;
    const double dvx = p->vx - planet->vx;
    const double dvy = p->vy - planet->vy;
    const double dvz = p->vz - planet->vz;
    const double r2 = dx*dx+dy*dy+dz*dz;
    const double d = sqrt (r2);
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    const double vcirc2 = mu/d;
    const double diff2 = v2 - vcirc2;
    const double a0 = -mu/( v2 - 2.*vcirc2 );

    const double j2 = (mag_p * mag_p * planet->r * planet->r * planet->r) / (3. * planet->m) * k2;
    const double lr = pow(2.*j2 * planet->r * planet->r * ap*ap*ap * pow((1 - ep*ep),(3./2.)) * planet->m / primary->m, (1./5.)); // star mass hard coded

    const double le_theta = theta_p - 0.5 * atan2(sin(2.*theta_p),cos(2.*theta_p) + (lr/a0)*(lr/a0)*(lr/a0)*(lr/a0)*(lr/a0)); // In the planet frame
    struct reb_vec3d beta_vec = reb_tools_spherical_to_xyz(1.,le_theta, phi_p);// In the planet frame
    reb_vec3d_irotate(&beta_vec, invrot); // rotates into xyz frame.
    const double v_dot_beta = dvz * beta_vec.x + dvy * beta_vec.y + dvz * beta_vec.z;
    //const double G = sim->G;
    struct reb_vec3d a = {0};
    if (tau_i < INFINITY){
        a.x = -2. * v_dot_beta * beta_vec.x / tau_i;
        a.y = -2. * v_dot_beta * beta_vec.y / tau_i;
        a.z = -2. * v_dot_beta * beta_vec.z / tau_i;
    }
    return a;
}

void rebx_disk_damping(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    const int N_active = sim->N_active;

    double k2 = 0.0;
    double tau_i = INFINITY;
    double mag_p = 0.0;
    double theta_p;
    double phi_p;
    int err=0;
    struct reb_rotation invrot;

    // For now just hard code in which planet we care about
    struct reb_particle* star = &particles[0];
    struct reb_particle* binary = &particles[1];
    const double* const k2_ptr = rebx_get_param(sim->extras, binary->ap, "k2");
    const double* const tau_i_ptr = rebx_get_param(sim->extras, binary->ap, "tau_inc_damping");
    const struct reb_vec3d* const Omega_ptr = rebx_get_param(sim->extras, binary->ap, "Omega");

    if (tau_i_ptr != NULL){
        tau_i = *tau_i_ptr;
    }
    if (k2_ptr != NULL){
        k2 = *k2_ptr;
    }
    if (Omega_ptr != NULL){
        struct reb_vec3d Omega = *Omega_ptr;
        struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *planet, *star, &err);

        struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, o.hvec);
        struct reb_rotation rot = reb_rotation_init_to_new_axes(o.hvec, line_of_nodes); // Invariant plane to planet
        if (isnan(rot.r)) {
          rot = reb_rotation_identity();
        }
        struct reb_vec3d Omega_rot = reb_vec3d_rotate(Omega, rot);
        reb_tools_xyz_to_spherical(Omega_rot, &mag_p, &theta_p, &phi_p); // In the frame of the planet, this motivates theta_p. Need this for calculating beta_test

        invrot = reb_rotation_inverse(rot);
    }

    // This is hard coded assuming all particles orbit the first planet
    for (int i = N_active; i < N; i++){
      struct reb_particle* test = &particles[i];
      struct reb_vec3d tot_force = rebx_calculate_disk_damping(sim, test, planet, star, tau_i, k2, mag_p, theta_p, phi_p, invrot);
      test->ax += tot_force.x;
      test->ay += tot_force.y;
      test->az += tot_force.z;
      //printf("Laplace Damping happening\n");
    }

}
