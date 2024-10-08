/** * @file solar_force.c
 * @brief   Add basic solar force
 * @author  Max Goldberg <max.goldberg@oca.eu>
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
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_solar_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double t, const double a_sun, const double P_sun, const double inc_sun, const double Omega_sun, const double m_sun, const int source_index){
    const struct reb_particle source = particles[source_index];
    const double G = sim->G;
    const double l_sun = 2*M_PI*t/P_sun;
    const double x0 = a_sun*cos(l_sun - Omega_sun);
    const double y0 = a_sun*sin(l_sun - Omega_sun);
    const double x = x0*cos(Omega_sun) - y0*sin(Omega_sun)*cos(inc_sun) + source.x;
    const double y = x0*sin(Omega_sun) + y0*cos(Omega_sun)*cos(inc_sun) + source.y;
    const double z = y0*sin(inc_sun) + source.z;
    const double r2_sun_primary = (x-source.x)*(x-source.x) + (y-source.y)*(y-source.y) + (z-source.z)*(z-source.z);
    const double r_sun_primary = sqrt(r2_sun_primary);
    const double fac_primary = -G*m_sun/r_sun_primary/r2_sun_primary;
    const double ax_on_primary = fac_primary*(source.x-x);
    const double ay_on_primary = fac_primary*(source.y-y);
    const double az_on_primary = fac_primary*(source.z-z);
    for (int i=0; i<N; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double r2_sun_moon = (x-p.x)*(x-p.x) + (y-p.y)*(y-p.y) + (z-p.z)*(z-p.z);
        const double r_sun_moon = sqrt(r2_sun_moon);
        const double fac_secondary = -G*m_sun/r_sun_moon/r2_sun_moon;
        const double ax_on_moon = fac_secondary*(p.x-x);
        const double ay_on_moon = fac_secondary*(p.y-y);
        const double az_on_moon = fac_secondary*(p.z-z);
        particles[i].ax += ax_on_moon - ax_on_primary;
        particles[i].ay += ay_on_moon - ay_on_primary;
        particles[i].az += az_on_moon - az_on_primary;
    }
}

void rebx_solar_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    const double* const inc = rebx_get_param(rebx, force->ap, "sf_inc");
    const double* const Omega = rebx_get_param(rebx, force->ap, "sf_Omega");
    const double* const a = rebx_get_param(rebx, force->ap, "sf_a");
    const double* const P_sun = rebx_get_param(rebx, force->ap, "sf_P");
    const double* const m = rebx_get_param(rebx, force->ap, "sf_m");
    if (inc != NULL && Omega != NULL && a != NULL && m != NULL){
        rebx_calculate_solar_force(sim, particles, N, sim->t, *a, *P_sun, *inc, *Omega, *m, 0); 
    }
}
