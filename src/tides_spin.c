/**
 * @file    tides_spin.c
 * @brief   Add self-consistent spin, tidal and dynamical equations of motion for bodies with structure
 * @author  Tiger Lu <tiger.lu@yale.edu>, Dan Tamayo <tamayo.daniel@gmail.com>, Hanno, Rein <hanno.rein@utoronto.ca>
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
 * $Tides$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Tiger Lu, Hanno Rein, D. Tamayo, Sam Hadden, Rosemary Mardling, Sarah Millholland, Gregory Laughlin
 * Implementation Paper    `Lu et al., 2023 <https://arxiv.org/abs/2303.00006>`_.
 * Based on                `Eggleton et al. 1998 <https://ui.adsabs.harvard.edu/abs/1998ApJ...499..853E/abstract>`_.
 * C Example               :ref:`c_example_tides_spin_pseudo_synchronization`, :ref:`c_example_tides_spin_migration_driven_obliquity_tides`, :ref:`c_example_tides_spin_kozai`.
 * Python Example          `SpinsIntro.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/SpinsIntro.ipynb>`_, `TidesSpinPseudoSynchronization.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesSpinPseudoSynchronization.ipynb>`_, `TidesSpinEarthMoon.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesSpinEarthMoon.ipynb>`_.
 * ======================= ===============================================
 *
 * This effect consistently tracks both the spin and orbital evolution of bodies under constant-time lag tides raised on both the primary and on the orbiting bodies.
 * In all cases, we need to set masses for all the particles that will feel these tidal forces. Particles with only mass are point particles.
 *
 * Particles are assumed to have structure (i.e - physical extent & distortion from spin) if the following parameters are set: physical radius particles[i].r, potential Love number of degree 2 k2 (Q/(1-Q) in Eggleton 1998), and the spin angular rotation frequency vector Omega.
 * If we wish to evolve a body's spin components, the fully dimensional moment of inertia I must be set as well. If this parameter is not set, the spin components will be stationary. Note that if the body is a test particle, this is assumed to be the specific moment of inertia.
 * Finally, if we wish to consider the effects of tides raised on a specific body, we must set the constant time lag tau as well.
 *
 * For spins that are synchronized with a circular orbit, the constant time lag can be related to the tidal quality factor Q as tau = 1/(2*n*tau), with n the orbital mean motion.
 * See Lu et. al (in review) and Eggleton et. al (1998) above for discussion.
 *
 *
 * **Effect Parameters**
 *
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * particles[i].r (float)       Yes         Physical radius (required for contribution from tides raised on the body).
 * k2 (float)                   Yes         Potential Love number of degree 2.
 * Omega (reb_vec3d)            Yes         Angular rotation frequency (Omega_x, Omega_y, Omega_z)
 * I (float)                    No          Moment of inertia (for test particles, assumed to be the specific MoI I/m)
 * tau (float)                  No          Constant time lag. If not set, defaults to 0
 * ============================ =========== ==================================================================
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

// eccentricity functions
double rebx_h1(double e){
  return 1. + (3./2.)*(e*e) + (1./8.)*(e*e*e*e) / pow((1. - e*e), (9./2.));
}

double rebx_h2(double e){
  return 1. + (9./2.)*(e*e) + (5./8.)*(e*e*e*e) / pow((1. - e*e), (9./2.));
}

double rebx_h3(double e){
  return ((1. + 3. * (e*e) + 3. * (e*e*e*e))/8.) / pow((1. - e*e), (9./2.));
}

double rebx_h4(double e){
  return (1. + (15. * (e*e)/2.) + (45. * (e*e*e*e)/8.) * (5. * (e*e*e*e*e*e)/16.)) / pow((1. - e*e), 6.);
}

double rebx_h5(double e){
  return (1. + (31. * (e*e)/2.) + (255. * (e*e*e*e)/8.) + (185. * (e*e*e*e*e*e)/16.) + (25. * (e*e*e*e*e*e*e*e)/64.)) / pow((1. - e*e), 15./2.);
}

double rebx_curlyL(double r){
  double A = -11.7039;
  double B = 12.8994;
  double C = -2.12471;

  double rrj = log10(r/4.676e-4);

  double exponent = A * rrj + B * rrj + C;
  double conversion = 1.12234e-11;

  return pow(10., exponent) * conversion;
}

double rebx_calculate_radius_inflation(struct reb_particle* source, struct reb_particle* target, const double G, const double k2, const double tau, const struct reb_vec3d Omega, const double dOmega, const double I){
  const double ms = source->m;
  const double Rs = source->r;
  const double mt = target->m;
  const double mtot = ms + mt;
  const double mu_ij = ms * mt / mtot; // have already checked for 0 and inf
  const double big_a = k2 * (Rs * Rs * Rs * Rs * Rs);

  double magp = sqrt(Omega.x*Omega.x + Omega.y*Omega.y + Omega.z*Omega.z);

  // distance vector FROM j TO i
  const double dx = source->x - target->x;
  const double dy = source->y - target->y;
  const double dz = source->z - target->z;
  struct reb_vec3d dvec = {.x=dx,.y=dy,.z=dz};
  const double d2 = dx * dx + dy * dy + dz * dz;
  const double dr = sqrt(d2);

  const double dvx = source->vx - target->vx;
  const double dvy = source->vy - target->vy;
  const double dvz = source->vz - target->vz;
  struct reb_vec3d vvec = {.x=dvx,.y=dvy,.z=dvz};
  const double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
  const double dv = sqrt(dv2);

  struct reb_vec3d tidal_force = {0};

  double sigma = 4 * tau * G / (3. * Rs*Rs*Rs*Rs*Rs * k2);
  if (sigma != 0.0){
    // Eggleton et. al 1998 tidal (equation 45)
    const double d_dot_vel = dx*dvx + dy*dvy + dz*dvz;

    // first vector
    const double vec1_x = 3. * d_dot_vel * dx;
    const double vec1_y = 3. * d_dot_vel * dy;
    const double vec1_z = 3. * d_dot_vel * dz;

    // h vector - EKH
    const double hx = dy * dvz - dz * dvy;
    const double hy = dz * dvx - dx * dvz;
    const double hz = dx * dvy - dy * dvx;

    // h - r^2 Omega
    const double comp_2_x = hx - d2 * Omega.x;
    const double comp_2_y = hy - d2 * Omega.y;
    const double comp_2_z = hz - d2 * Omega.z;

    // second vector
    const double vec2_x = comp_2_y * dz - comp_2_z * dy;
    const double vec2_y = comp_2_z * dx - comp_2_x * dz;
    const double vec2_z = comp_2_x * dy - comp_2_y * dx;
    const double prefactor = (-9. * sigma * mt * mt * big_a * big_a) / (2. * mu_ij * (d2 * d2 * d2 * d2 * d2));

    tidal_force.x = (prefactor * (vec1_x + vec2_x));
    tidal_force.y = (prefactor * (vec1_y + vec2_y));
    tidal_force.z = (prefactor * (vec1_z + vec2_z));
  }

  double alpha = 0.261;// Hard coded for now
  double denominator = 0.5 * G * mt * mt / (Rs * Rs) + alpha * mt * Rs * magp * magp;

  // Curly L
  double curlyL = -1. * rebx_curlyL(Rs);

  // Tidal force
  struct reb_vec3d v3_1 = reb_vec3d_mul(tidal_force, mu_ij);
  struct reb_vec3d v3_22 = reb_vec3d_mul(reb_vec3d_cross(Omega, dvec), -1.);
  struct reb_vec3d v3_2 = reb_vec3d_add(vvec, v3_22);
  const double t3 = reb_vec3d_dot(v3_1, v3_1);

  double retval = 0.5 * (curlyL + t3) / denominator;
  return retval;
}

struct reb_vec3d rebx_calculate_spin_orbit_accelerations(struct reb_particle* source, struct reb_particle* target, const double G, const double k2, const double sigma, const struct reb_vec3d Omega){
  // All quantities associated with SOURCE
  // This is the quadrupole potential/tides raised on the SOURCE
  const double ms = source->m;
  const double Rs = source->r;
  const double mt = target->m;
  const double mtot = ms + mt;
  const double mu_ij = ms * mt / mtot; // have already checked for 0 and inf
  const double big_a = k2 * (Rs * Rs * Rs * Rs * Rs);

  // distance vector FROM j TO i
  const double dx = source->x - target->x;
  const double dy = source->y - target->y;
  const double dz = source->z - target->z;
  const double d2 = dx * dx + dy * dy + dz * dz;
  const double dr = sqrt(d2);

  // Velocity vector: i to j
  const double dvx = source->vx - target->vx;
  const double dvy = source->vy - target->vy;
  const double dvz = source->vz - target->vz;
  //const double vel2 = dvx * dvx + dvy * dvy + dvz * dvz;
  //const double vr = sqrt(vel2);

  struct reb_vec3d tot_force = {0};

  if (k2 != 0.0){
    // Eggleton et. al 1998 quadrupole (equation 33)
    const double quad_prefactor = mt * big_a / mu_ij;
    const double omega_dot_d = Omega.x * dx + Omega.y * dy + Omega.z * dz;
    const double omega_squared = Omega.x * Omega.x + Omega.y * Omega.y + Omega.z * Omega.z;

    const double t1 = 5. * omega_dot_d * omega_dot_d / (2. * (dr * dr * dr * dr * dr * dr * dr));
    const double t2 = omega_squared / (2. * (dr * dr * dr * dr * dr));
    const double t3 = omega_dot_d / (dr * dr * dr * dr * dr);
    const double t4 = 6. * G * mt / (dr * dr * dr * dr * dr * dr * dr * dr);

    tot_force.x = (quad_prefactor * ((t1 - t2 - t4) * dx - (t3 * Omega.x)));
    tot_force.y = (quad_prefactor * ((t1 - t2 - t4) * dy - (t3 * Omega.y)));
    tot_force.z = (quad_prefactor * ((t1 - t2 - t4) * dz - (t3 * Omega.z)));

    if (sigma != 0.0){
      // Eggleton et. al 1998 tidal (equation 45)
      const double d_dot_vel = dx*dvx + dy*dvy + dz*dvz;

      // first vector
      const double vec1_x = 3. * d_dot_vel * dx;
      const double vec1_y = 3. * d_dot_vel * dy;
      const double vec1_z = 3. * d_dot_vel * dz;

      // h vector - EKH
      const double hx = dy * dvz - dz * dvy;
      const double hy = dz * dvx - dx * dvz;
      const double hz = dx * dvy - dy * dvx;

      // h - r^2 Omega
      const double comp_2_x = hx - d2 * Omega.x;
      const double comp_2_y = hy - d2 * Omega.y;
      const double comp_2_z = hz - d2 * Omega.z;

      // second vector
      const double vec2_x = comp_2_y * dz - comp_2_z * dy;
      const double vec2_y = comp_2_z * dx - comp_2_x * dz;
      const double vec2_z = comp_2_x * dy - comp_2_y * dx;

      const double prefactor = (-9. * sigma * mt * mt * big_a * big_a) / (2. * mu_ij * (d2 * d2 * d2 * d2 * d2));

      tot_force.x += (prefactor * (vec1_x + vec2_x));
      tot_force.y += (prefactor * (vec1_y + vec2_y));
      tot_force.z += (prefactor * (vec1_z + vec2_z));
    }
  }

  return tot_force;
}

static void rebx_spin_orbit_accelerations(struct reb_particle* source, struct reb_particle* target, const double G, const double k2, const double sigma, const struct reb_vec3d Omega){

    // Input params all associated with source
    const double ms = source->m;
    const double mt = target->m;
    const double mtot = ms + mt;

    // check if ODE is set here
    struct reb_vec3d tot_force = rebx_calculate_spin_orbit_accelerations(source, target, G, k2, sigma, Omega);

    target->ax -= ((ms / mtot) * tot_force.x);
    target->ay -= ((ms / mtot) * tot_force.y);
    target->az -= ((ms / mtot) * tot_force.z);

    source->ax += ((mt / mtot) * tot_force.x);
    source->ay += ((mt / mtot) * tot_force.y);
    source->az += ((mt / mtot) * tot_force.z);

}

static void rebx_spin_derivatives(struct reb_ode* const ode, double* const yDot, const double* const y, const double t){
    struct reb_simulation* sim = ode->ref;
    struct rebx_extras* const rebx = sim->extras;
    unsigned int Nspins = 0;
    const int N_real = sim->N - sim->N_var;
    for (int i=0; i<N_real; i++){
        struct reb_particle* pi = &sim->particles[i]; // target particle
        const double* k2 = rebx_get_param(rebx, pi->ap, "k2"); // This is slow
        const double* tau = rebx_get_param(rebx, pi->ap, "tau");
        const double* I = rebx_get_param(rebx, pi->ap, "I");

        // Particle MUST have k2 and moment of inertia to feel effects
        if (k2 != NULL && I != NULL){

          // Tidal dissipation off by default. Check for non-zero tau here.
          double sigma_in = 0.0;
          if (tau != NULL){
            sigma_in = 4 * (*tau) * sim->G / (3. * pi->r * pi->r * pi->r * pi->r * pi->r * (*k2));
          }
	       // Set initial spin accelerations to 0
          yDot[4*Nspins] = 0;
          yDot[4*Nspins + 1] = 0;
          yDot[4*Nspins + 2] = 0;
          yDot[4*Nspins + 3] = 0;

          const struct reb_vec3d Omega = {.x=y[4*Nspins], .y=y[4*Nspins+1], .z=y[4*Nspins+2]};
          for (int j=0; j<N_real; j++){
            if (i != j){
                struct reb_particle* pj = &sim->particles[j];

                const double mi = pi->m;
                const double mj = pj->m;
                double mu_ij;

            		if (mj == 0){
            		   continue;
            		}

                double I_specific;
                if (mi == 0){ // If test particle, assume I = specific moment of inertia
                    I_specific = *I;
                }
                else{
                  mu_ij = (mi * mj) / (mi + mj);
                  I_specific = *I / mu_ij;
                }

		// di - dj
                const double dx = pi->x - pj->x;
                const double dy = pi->y - pj->y;
                const double dz = pi->z - pj->z;

		            struct reb_vec3d tf = rebx_calculate_spin_orbit_accelerations(pi, pj, sim->G, *k2, sigma_in, Omega);
                // Eggleton et. al 1998 spin EoM (equation 36)
                yDot[4*Nspins] += ((dy * tf.z - dz * tf.y) / (-I_specific));
                yDot[4*Nspins + 1] += ((dz * tf.x - dx * tf.z) / (-I_specific));
                yDot[4*Nspins + 2] += ((dx * tf.y - dy * tf.x) / (-I_specific));
            }
          }

          // Radius inflation
          if (i != 0){
            double dOmega = sqrt(yDot[4*Nspins]*yDot[4*Nspins]+yDot[4*Nspins+1]*yDot[4*Nspins+1]+yDot[4*Nspins+2]*yDot[4*Nspins+2]);
            yDot[4*Nspins + 3] = rebx_calculate_radius_inflation(pi, &sim->particles[0], sim->G, *k2, *tau, Omega, dOmega, *I);//
            //printf("%d %e\n", i, yDot[4*Nspins + 3]);
            //exit(1);
          }
          Nspins += 1;
      }
    }
    if (ode->length != Nspins*4){
        reb_simulation_error(sim, "rebx_spin ODE is not of the expected length.\n");
        exit(1);
    }
    //exit(1);

}

static void rebx_spin_sync_pre(struct reb_ode* const ode, const double* const y0){
    struct reb_simulation* sim = ode->ref;
    struct rebx_extras* const rebx = sim->extras;
    unsigned int Nspins = 0;
    const int N_real = sim->N - sim->N_var;
    for (int i=0; i<N_real; i++){
        struct reb_particle* p = &sim->particles[i];
        double* I = rebx_get_param(rebx, p->ap, "I");
        struct reb_vec3d* Omega = rebx_get_param(rebx, p->ap, "Omega");
        if (I != NULL && Omega != NULL){
            const struct reb_vec3d* Omega = rebx_get_param(rebx, p->ap, "Omega");
            ode->y[4*Nspins] = Omega->x;
            ode->y[4*Nspins+1] = Omega->y;
            ode->y[4*Nspins+2] = Omega->z;

            // radius
            ode->y[4*Nspins+3] = p->r;
            Nspins += 1;
        }
    }

    if (ode->length != Nspins*4){
        reb_simulation_error(sim, "rebx_spin ODE is not of the expected length.\n");
        exit(1);
    }
}

static void rebx_spin_sync_post(struct reb_ode* const ode, const double* const y0){
    struct reb_simulation* sim = ode->ref;
    struct rebx_extras* const rebx = sim->extras;
    unsigned int Nspins = 0;
    const int N_real = sim->N - sim->N_var;
    for (int i=0; i<N_real; i++){
        struct reb_particle* p = &sim->particles[i];
        double* I = rebx_get_param(rebx, p->ap, "I");
        struct reb_vec3d* Omega = rebx_get_param(rebx, p->ap, "Omega");
        if (I != NULL && Omega != NULL){
            rebx_set_param_vec3d(rebx, (struct rebx_node**)&p->ap, "Omega", (struct reb_vec3d){.x=y0[4*Nspins], .y=y0[4*Nspins+1], .z=y0[4*Nspins+2]});

            // radius
            p->r = ode->y[4*Nspins+3];
            Nspins += 1;
        }
    }
    if (ode->length != Nspins*4){
        reb_simulation_error(sim, "rebx_spin ODE is not of the expected length.\n");
        exit(0);
    }
}

void rebx_spin_initialize_ode(struct rebx_extras* const rebx, struct rebx_force* const effect){
    struct reb_simulation* sim = rebx->sim;
    unsigned int Nspins = 0;
    const int N_real = sim->N - sim->N_var;
    for (int i=0; i<N_real; i++){
        struct reb_particle* p = &sim->particles[i];
        // Only track spin if particle has moment of inertia and valid spin axis set
        double* I = rebx_get_param(rebx, p->ap, "I");
        struct reb_vec3d* Omega = rebx_get_param(rebx, p->ap, "Omega");
        if (I != NULL && Omega != NULL){
            Nspins += 1;
        }
    }

    if (Nspins > 0){
        struct reb_ode* spin_ode = reb_ode_create(sim, Nspins*4);
        spin_ode->ref = sim;
        spin_ode->derivatives = rebx_spin_derivatives;
        spin_ode->pre_timestep = rebx_spin_sync_pre;
        spin_ode->post_timestep = rebx_spin_sync_post;
        rebx_set_param_pointer(rebx, &effect->ap, "ode", spin_ode);
    }
}

void rebx_tides_spin(struct reb_simulation* const sim, struct rebx_force* const effect, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    const double G = sim->G;

    // check if ODE is initialized
    struct reb_ode** ode = sim->odes;
    if (ode == NULL){
      reb_simulation_warning(sim, "Spin axes are not being evolved. Call rebx_spin_initialize_ode to evolve\n");
    }

    for (int i=0; i<N; i++){
        struct reb_particle* source = &particles[i];
        // Particle must have a k2 set, otherwise we treat this body as a point particle
        const double* k2 = rebx_get_param(rebx, source->ap, "k2");
        const double* tau = rebx_get_param(rebx, source->ap, "tau");
        const struct reb_vec3d* Omega = rebx_get_param(rebx, source->ap, "Omega");

        // Particle needs all three spin components and k2 to feel additional forces
        if (Omega != NULL && k2 != NULL){
          // Tidal dissipation off by default. Check for non-zero tau here.
          double sigma_in = 0.0;
          if (tau != NULL){
            sigma_in = 4 * (*tau) * sim->G / (3. * source->r * source->r * source->r * source->r * source->r * (*k2));
          }

          for (int j=0; j<N; j++){
              if (i==j){
                  continue;
              }
              struct reb_particle* target = &particles[j]; // j raises tides on i
              if (source->m == 0 || target->m == 0){
                  continue;
              }

              rebx_spin_orbit_accelerations(source, target, G, *k2, sigma_in, *Omega);
          }
      }
    }
}

// Calculate potential of conservative piece of interaction between a point mass target and a source with a tidally and rotationally induced quadrupole
// Equation 31 in Eggleton et. al (1998)
static double rebx_calculate_spin_potential(struct reb_particle* source, struct reb_particle* target, const double G, const double k2, const struct reb_vec3d Omega){
    const double Rs = source->r;
    const double mt = target->m;

    const double big_a = k2 * (Rs * Rs * Rs * Rs * Rs);

    // distance vector FROM j TO i
    const double dx = source->x - target->x;
    const double dy = source->y - target->y;
    const double dz = source->z - target->z;
    const double d2 = dx * dx + dy * dy + dz * dz;
    const double dr = sqrt(d2);

    const double omega_dot_d = Omega.x * dx + Omega.y * dy + Omega.z * dz;
    const double omega_squared = Omega.x * Omega.x + Omega.y * Omega.y + Omega.z * Omega.z;

    const double t1 = -mt * big_a * omega_dot_d * omega_dot_d / (2 * d2 * d2 * dr);
    const double t2 = mt * big_a * omega_squared / (6 * d2 * dr);
    const double t3 = G * mt * mt * big_a / (d2 * d2 * d2);

    return -(t1 + t2 + t3);
}

double rebx_tides_spin_energy(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    struct reb_simulation* const sim = rebx->sim;
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    const double G = sim->G;
    double E=0.;

    for (int i=0; i<N_real; i++){
        struct reb_particle* source = &particles[i];
        // Particle must have a k2, radius and mass set, otherwise we treat this body as a point particle
        const double* k2 = rebx_get_param(rebx, source->ap, "k2");
        const struct reb_vec3d* Omegaptr = rebx_get_param(rebx, source->ap, "Omega");
        if (k2 == NULL || source->m == 0 || source->r == 0){
            continue;
        }
        struct reb_vec3d Omega = {0};
        if (Omegaptr != NULL){
            Omega = *Omegaptr;
        }
        double* I = rebx_get_param(rebx, source->ap, "I");
        if (I != NULL){
            const double omega_squared = Omega.x * Omega.x + Omega.y * Omega.y + Omega.z * Omega.z;
            E += 0.5 * (*I) * omega_squared;
        }
        for (int j=0; j<N_real; j++){
            if (i==j){
                continue;
            }
            struct reb_particle* target = &particles[j]; // planet raising the tides on the star
            if (target->m > 0){
                E += rebx_calculate_spin_potential(source, target, G, *k2, Omega);
            }
        }
    }

    return E;
}
