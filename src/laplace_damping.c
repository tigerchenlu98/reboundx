/**
 * @file    type_I_migration.c
 * @brief   Type I migration
 * @author  Kaltrina Kajtazi <1kaltrinakajtazi@gmail.com>, Gabriele Pichierri <gabrielepichierri@gmail.com>
 *
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Orbit Modifications$       // Effect category
 *
 * ======================= ===============================================
 * Authors                 Kajtazi, Kaltrina and D. Petit, C. Antoine
 * Implementation Paper    `Kajtazi et al 2022 <https://ui.adsabs.harvard.edu/abs/2022arXiv221106181K/abstract>`_.
 * Based on                `Cresswell & Nelson 2008 <https://ui.adsabs.harvard.edu/abs/2008A%26A...482..677C/abstract>`_, and `Pichierri et al 2018 <https://ui.adsabs.harvard.edu/abs/2018CeMDA.130...54P/abstract>`_.
 * C example               :ref:`c_example_type_I_migration`
 * Python example          `TypeIMigration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TypeIMigration.ipynb>`_.
 * ======================= ===============================================
 *
 * This applies Type I migration, damping eccentricity, angular momentum and inclination.
 * The base of the code is the same as the modified orbital forces one written by D. Tamayo, H. Rein.
 * It also allows for parameters describing an inner disc edge, modeled using the implementation in inner_disk_edge.c.
 * Note that this code is not machine independent since power laws were not possible to avoid all together.
 *
 * **Effect Parameters**
 *
 * ===================================== =========== ==================================================================================================================
 * Field (C type)                        Required    Description
 * ===================================== =========== ==================================================================================================================
 * ide_position (double)                 No          The position of the inner disk edge in code units
 * ide_width (double)                    No          The disk edge width (planet will stop within ide_width of ide_position)
 * tIm_surface_density_1 (double)        Yes         Disk surface density at one code unit from the star; used to find the surface density at any distance from the star
 * tIm_scale_height_1 (double)           Yes         The scale height at one code unit from the star; used to find the aspect ratio at any distance from the star
 * tIm_surface_density_exponent (double) Yes         Exponent of disk surface density, indicative of the surface density profile of the disk
 * tIm_flaring_index (double)            Yes         The flaring index; 1 means disk is irradiated by only the stellar flux
 * ===================================== =========== ==================================================================================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static struct reb_vec3d rebx_calculate_laplace_damping(struct reb_simulation* const sim, struct reb_particle* p, struct reb_particle* planet, struct reb_particle* primary, double tau_i, double k2, double mag_p, double theta_p, double phi_p, struct reb_rotation invrot){

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
    const double lr = pow(j2 * planet->r * planet->r * ap*ap*ap * pow((1 - ep*ep),(3./2.)) * planet->m / primary->m, (1./5.)); // star mass hard coded

    const double le_theta = theta_p - 0.5 * atan2(sin(2.*theta_p),cos(2.*theta_p) + (lr/a0)*(lr/a0)*(lr/a0)*(lr/a0)*(lr/a0)); // In the planet frame
    struct reb_vec3d beta_vec = reb_tools_spherical_to_xyz(1.,le_theta, phi_p);// In the planet frame
    reb_vec3d_irotate(&beta_vec, invrot); // rotates into xya frame.
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

void rebx_laplace_damping(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
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
    struct reb_particle* planet = &particles[1];
    const double* const k2_ptr = rebx_get_param(sim->extras, planet->ap, "k2");
    const double* const tau_i_ptr = rebx_get_param(sim->extras, planet->ap, "tau_i");
    const struct reb_vec3d* const Omega_ptr = rebx_get_param(sim->extras, planet->ap, "Omega");

    if (tau_i_ptr != NULL){
        tau_i = *tau_i_ptr;
    }
    if (k2_ptr != NULL){
        k2 = *k2_ptr;
    }
    if (Omega_ptr != NULL){
        struct reb_vec3d Omega = *Omega_ptr;
        struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *planet, *star, &err);
        struct reb_vec3d h = o.hvec;
        struct reb_vec3d xyz = {0., 0., 1.};

        struct reb_rotation rot = reb_rotation_init_from_to(xyz, h); // rotates from xyz to planet frame
        invrot = reb_rotation_inverse(rot); // rotates from planet frame to xyz
        struct reb_vec3d Omega_rot = reb_vec3d_rotate(Omega, rot);
        reb_tools_xyz_to_spherical(Omega_rot, &mag_p, &theta_p, &phi_p); // In the frame of the planet, this motivates theta_p
    }

    // This is hard coded assuming all particles orbit the first planet
    for (int i = N_active; i < N; i++){
      struct reb_particle* test = &particles[i];
      struct reb_vec3d tot_force = rebx_calculate_laplace_damping(sim, test, planet, star, tau_i, k2, mag_p, theta_p, phi_p, invrot);
      test->ax += tot_force.x;
      test->ay += tot_force.y;
      test->az += tot_force.z;
    }

}
