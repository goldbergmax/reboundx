/**
 * @file    turbulent_force.c
 * @brief   Add turbulent forces
 * @author  Max Goldberg <max.goldberg@oca.eu>
 * 
 * @section     LICENSE
 * Copyright (c) 2022 Hanno Rein, Dan Tamayo
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
 * $Turbulent Forces$     
 *
 * ======================= ===============================================
 * Authors                 H. Rein
 * Based on                `Baruteau and Lin 2010 <https://ui.adsabs.harvard.edu/abs/2010ApJ...709..759B>`_.
 * C Example               :ref:`c_example_stochastic_forces`
 * Python Example          `StochasticForces.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/StochasticForces.ipynb>`_, `StochasticForcesCartesian.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/StochasticForcesCartesian.ipynb>`_,
 * ======================= ===============================================
 * 
 * This manages a set of turbulent modes and applies turbulent forces to particles in the simulation.  
 * 
 * **Effect Parameters**
 * 
 * The particle with index 0 cannot experience turbulent forces.
 * 
 * ==================================== =========== ==================================================================================
 * Field (C type)                       Required    Description
 * ==================================== =========== ==================================================================================
 * turb_gamma (double)                  Yes         Strength of turbulence 
 * turb_Gamma (double)                  Yes         Response factor of disk
 * turb_h0 (double)                     Yes         Disk scale height at 1 AU
 * turb_flaring_index (double)          Yes         Flaring index of disk
 * turb_modes (struct rebx_turb_modes*) Yes         Pointer to struct containing turbulent modes
 * ==================================== =========== ==================================================================================
 *
 * **Particle Parameters**
 * 
 * None
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"


static void rebx_random_normal2(struct reb_simulation* r, double* n0, double* n1){
	double v1,v2,rsq=1.;
	while(rsq>=1. || rsq<1.0e-12){
		v1=2.*((double)rand_r(&(r->rand_seed)))/((double)(RAND_MAX))-1.0;
		v2=2.*((double)rand_r(&(r->rand_seed)))/((double)(RAND_MAX))-1.0;
		rsq=v1*v1+v2*v2;
	}
	*n0 = v1*sqrt(-2.*log(rsq)/rsq);
	*n1 = v2*sqrt(-2.*log(rsq)/rsq);
}

void rebx_update_modes(struct reb_simulation* const sim, struct rebx_force* const force){
    struct rebx_extras* const rebx = sim->extras;
    struct rebx_turb_modes* modes = rebx_get_param(rebx, force->ap, "turb_modes");
    double *h0 = rebx_get_param(rebx, force->ap, "turb_h0");
    double *flaring_index = rebx_get_param(rebx, force->ap, "turb_flaring_index");
    double *inner_edge = rebx_get_param(rebx, force->ap, "turb_inner_edge");
    double *outer_edge = rebx_get_param(rebx, force->ap, "turb_outer_edge");
    int *max_m_ptr = rebx_get_param(rebx, force->ap, "turb_max_m");
    if (modes == NULL){
        reb_simulation_error(sim, "No modes found in rebx_turbulent_forces.\n");
        return;
    }
    if (h0 == NULL || flaring_index == NULL){
        reb_simulation_error(sim, "h0 or flaring_index not set in rebx_turbulent_forces.\n");
        return;
    }
    if (inner_edge == NULL || outer_edge == NULL){
        reb_simulation_error(sim, "inner_edge or outer_edge not set in rebx_turbulent_forces.\n");
        return;
    }
    int max_m = 96;
    if (max_m_ptr != NULL){
        max_m = *max_m_ptr;
        if (max_m < 0){
            reb_simulation_error(sim, "max_m must be positive in rebx_turbulent_forces.\n");
            return;
        }
    }
    for (int i=0; i<modes->nmodes; i++){
        if (sim->t >= modes->t0k[i] + modes->Deltatk[i]){
            modes->xi[i] = reb_random_normal(sim, 1.0);
            modes->m[i] = (int)pow(10, reb_random_uniform(sim, 0.0, log10(max_m)));
            modes->rc[i] = reb_random_uniform(sim, *inner_edge, *outer_edge);
            modes->phik[i] = reb_random_uniform(sim, 0.0, 2*M_PI);
            modes->sigmak[i] = M_PI*modes->rc[i]/(4*((double)modes->m[i]));
            modes->Omegak[i] = sqrt(sim->G*sim->particles[0].m/(modes->rc[i]*modes->rc[i]*modes->rc[i]));
            const double cs = (*h0) * pow(modes->rc[i], (*flaring_index)) * modes->rc[i] * modes->Omegak[i];
            modes->Deltatk[i] = 2*M_PI/(((double)modes->m[i])*cs);
            modes->t0k[i] = sim->t;
        }
    }
}

void rebx_turbulent_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    struct reb_particle com = particles[0];
    rebx_update_modes(sim, force);
    struct rebx_turb_modes* modes = rebx_get_param(rebx, force->ap, "turb_modes");
    double *gamma = rebx_get_param(rebx, force->ap, "turb_gamma");
    double *Gamma = rebx_get_param(rebx, force->ap, "turb_Gamma");
    if (gamma == NULL || Gamma == NULL){
        reb_simulation_error(sim, "gamma or Gamma not set in rebx_turbulent_forces.\n");
        return;
    }
    for (int i=0; i<N; i++){
        if (modes == NULL){
            reb_simulation_error(sim, "No modes found in rebx_turbulent_forces.\n");
            return;
        }
        if (i>0 && modes != NULL){
            const struct reb_particle p = particles[i];

            int err=0;
            struct reb_orbit o = reb_orbit_from_particle_err(sim->G, particles[i], com, &err);
            if (err){
                reb_simulation_error(sim, "An error occured during the orbit calculation in rebx_stochastic_forces.\n");
                return;
            }

            double turbulent_force_r = 0.0;
            double turbulent_force_phi = 0.0;

            const double dx = p.x - com.x; 
            const double dy = p.y - com.y;
            const double dz = p.z - com.z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);
            const double phi = atan2(dy, dx);
            
            const double force_prefac = (*gamma)*(*Gamma)*dr*o.n*o.n;
            for (int j=0; j<modes->nmodes; j++){
                const double expfac = exp(-(dr - modes->rc[j])*(dr - modes->rc[j])/(modes->sigmak[j]*modes->sigmak[j]));
                const double arg1 = ((double)modes->m[j])*phi - modes->phik[j] - modes->Omegak[j]*(sim->t - modes->t0k[j]); // check with morby about m multiplying everything
                const double arg2 = M_PI*(sim->t - modes->t0k[j])/modes->Deltatk[j];
                const double r_prefac = (1+2*dr*(dr - modes->rc[j])/(modes->sigmak[j]*modes->sigmak[j]));
                turbulent_force_r += r_prefac*modes->xi[j]*expfac*cos(arg1)*sin(arg2);
                turbulent_force_phi += ((double)modes->m[j])*modes->xi[j]*expfac*sin(arg1)*sin(arg2);
            }

            particles[i].ax += force_prefac*(turbulent_force_r*dx/dr - turbulent_force_phi*dy/dr);
            particles[i].ay += force_prefac*(turbulent_force_r*dy/dr + turbulent_force_phi*dx/dr);


		    com = reb_particle_com_of_pair(com, p);
        }
    }
}
