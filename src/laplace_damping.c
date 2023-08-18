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

static struct reb_vec3d rebx_calculate_laplace_damping(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){

    /* Default values for the parameters in case the user forgets to define them when using this code */
    double dedge = 0.0;
    double hedge = 0.0;
    double invtau_mig = 0.0;
    double tau_i = INFINITY;
    double lr = 0.0;
    double mag_p;
    double theta_p;
    double phi_p;

    /* Parameters that should be changed/set in Python notebook or in C outside of this */
    const double* const dedge_ptr = rebx_get_param(sim->extras, force->ap, "ide_position");
    const double* const hedge_ptr = rebx_get_param(sim->extras, force->ap, "ide_width");

    const double* const lr_ptr = rebx_get_param(sim->extras, source->ap, "lr");
    const double* const tau_i_ptr = rebx_get_param(sim->extras, p->ap, "tau_i");
    const struct reb_vec3d* Omega = rebx_get_param(sim->extras, source->ap, "Omega");

    if (dedge_ptr != NULL){
        dedge = *dedge_ptr;
    }
    if (hedge_ptr != NULL){
        hedge = *hedge_ptr;
    }
    if (tau_i_ptr != NULL){
        tau_i = *tau_i_ptr;
    }
    if (lr_ptr != NULL){
        lr = *lr_ptr;
    }
    if (Omega != NULL){
        reb_tools_xyz_to_spherical(*Omega, &mag_p, &theta_p, &phi_p);
    }

    /* Accessing the calculated semi-major axis, eccentricity and inclination for each integration step, via modify_orbits_direct where they are calculated and returned*/
    int err=0;
    struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *source, &err);
    const double a0 = o.a;

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x-source->x;
    const double dy = p->y-source->y;
    const double dz = p->z-source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    const double le_theta = theta_p - 0.5 * atan2(sin(2.*theta_p),cos(2.*theta_p) + (lr/a0)*(lr/a0)*(lr/a0)*(lr/a0)*(lr/a0));
    struct reb_vec3d beta_vec = reb_tools_spherical_to_xyz(1.,le_theta, phi_p);
    const double v_dot_beta = dvz * beta_vec.x + dvy * beta_vec.y + dvz * beta_vec.z;

    /* Calculating the aspect ratio evaluated at the position of the planet, r and defining other variables */

    const double G = sim->G;
    struct reb_vec3d a = {0};
    if (tau_i < INFINITY){
        a.x = -2. * v_dot_beta * beta_vec.x / tau_i;
        a.y = -2. * v_dot_beta * beta_vec.y / tau_i;
        a.z = -2. * v_dot_beta * beta_vec.z / tau_i;
       //invtau_mig = rebx_calculate_planet_trap(a0, dedge, hedge) / tau_i;
    }

    //if (invtau_mig != 0.0){
    //    a.x = -dvx*(invtau_mig);
    //    a.y = -dvy*(invtau_mig);
    //    a.z = -dvz*(invtau_mig);
    //}


    return a;
}

void rebx_laplace_damping(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    //int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_PARTICLE; // Default
    //if (ptr != NULL){
    //    coordinates = *ptr;
    //}
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_laplace_damping, particles, N);
}
