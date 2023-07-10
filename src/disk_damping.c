/**
 * @file    modify_orbits_forces.c
 * @brief   Update orbital elements with prescribed timescales using forces.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
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
 * $Orbit Modifications$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo, H. Rein
 * Implementation Paper    `Kostov et al., 2016 <https://ui.adsabs.harvard.edu/abs/2016ApJ...832..183K/abstract>`_.
 * Based on                `Papaloizou & Larwood 2000 <http://labs.adsabs.harvard.edu/adsabs/abs/2000MNRAS.315..823P/>`_.
 * C Example               :ref:`c_example_modify_orbits`
 * Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_
 *                         `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
 * ======================= ===============================================
 *
 * This applies physical forces that orbit-average to give exponential growth/decay of the semimajor axis, eccentricity and inclination.
 * The eccentricity damping keeps the angular momentum constant (corresponding to `p=1` in modify_orbits_direct), which means that eccentricity damping will induce some semimajor axis evolution.
 * Additionally, eccentricity/inclination damping will induce pericenter/nodal precession.
 * Both these effects are physical, and the method is more robust for strongly perturbed systems.
 *
 * **Effect Parameters**
 *
 * If coordinates not, defaults to using Jacobi coordinates.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
 *                                          See the examples for usage.
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * One can pick and choose which particles have which parameters set.
 * For each particle, any unset parameter is ignored.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * tau_a (double)               No          Semimajor axis exponential growth/damping timescale
 * tau_e (double)               No          Eccentricity exponential growth/damping timescale
 * tau_inc (double)             No          Inclination axis exponential growth/damping timescale
 * ============================ =========== ==================================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static double rebx_acos2(double num, double denom, double disambiguator){ // Tlu hack put this somewhere else
	double val;
	double cosine = num/denom;
	if(cosine > -1. && cosine < 1.){
		val = acos(cosine);
		if(disambiguator < 0.){
			val = - val;
		}
	}
	else{
		val = (cosine <= -1.) ? M_PI : 0.;
	}
	return val;
}

static struct reb_particle rebx_calculate_disk_damping(struct reb_simulation* const sim, struct rebx_operator* const operator, struct reb_particle* p, struct reb_particle* primary, const double dt){
    int err=0;
    struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *primary, &err);
    if(err){        // mass of primary was 0 or p = primary.  Return same particle without doing anything.
        return *p;
    }
    double Cm = 0.0;
    double Ca = 0.0;
    double tc = INFINITY;

    const double* const Cm_ptr = rebx_get_param(sim->extras, p->ap, "dd_Cm");
    const double* const Ca_ptr = rebx_get_param(sim->extras, p->ap, "dd_Ca");
    const double* const tc_ptr = rebx_get_param(sim->extras, p->ap, "dd_tc");

    if(Cm_ptr != NULL){
        Cm = *Cm_ptr;
    }
    if(Ca_ptr != NULL){
        Ca = *Ca_ptr;
    }
    if(tc_ptr != NULL){
        tc = *tc_ptr;
    }

    if (Cm > 0.0 && tc < INFINITY){
        double de_dt = (-Cm * o.e / (2. * tc)) * (1. + (Ca / (o.e*o.e + o.inc*o.inc)));
        o.e += dt * de_dt;

        double di_dt = (-Cm * o.inc / (2. * tc)) * (1. + (Ca / (o.e*o.e + o.inc*o.inc)));
        o.inc += dt * di_dt;
        /*const double dvx = p->vx - primary->vx;
        const double dvy = p->vy - primary->vy;
        const double dvz = p->vz - primary->vz;
        const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
        const double dx = p->x-primary->x;
        const double dy = p->y-primary->y;
        const double dz = p->z-primary->z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double dr = sqrt(r2);

        const double mu = sim->G*(p->m+primary->m);
        const double muinv = 1./mu;
        const double vcircsquared = mu / dr;
        const double vdiffsquared = v2 - vcircsquared;

        const double vdotr = dx*dvx + dy*dvy + dz*dvz;

        double ex = muinv*( vdiffsquared*dx - vdotr*dvx );
        double ey = muinv*( vdiffsquared*dy - vdotr*dvy );
        double ez = muinv*( vdiffsquared*dz - vdotr*dvz );
        double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity

        double hx = (dy*dvz - dz*dvy); 					// specific angular momentum vector
        double hy = (dz*dvx - dx*dvz);
        double hz = (dx*dvy - dy*dvx);
        double h = sqrt( hx*hx + hy*hy + hz*hz );		// abs value of angular momentum
        double inc = rebx_acos2(hz, h, 1.); // inclination
        */
    }
    return reb_tools_orbit_to_particle(sim->G, *primary, p->m, o.a, o.e, o.inc, o.Omega, o.omega, o.f);
}

void rebx_disk_damping(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    const int* const ptr = rebx_get_param(sim->extras, operator->ap, "coordinates");
   	enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI;
  	if (ptr != NULL){
  		coordinates = *ptr;
  	}
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_tools_com_ptm(sim, operator, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_disk_damping, dt);
}
