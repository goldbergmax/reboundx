/**
 * @file    full_type_I.c
 * @brief   Type I migration 
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

// wave damping timescale from Tanaka & Ward 2004
const double rebx_calculate_wave_timescale(const double G, const double sd, const double r, const double ms, const double mp, const double sma, const double h2){    
    const double t_wave = (sqrt(ms*ms*ms)*h2*h2)/(mp*sd*sqrt(sma*G));
    return t_wave;
}

// Lindblad torque from PBK11 Eq. 3 and reduction factors from Coleman & Nelson 2014
const double rebx_calculate_lindblad_torque(const double sd_ind, const double temp_ind, const double adi_ind, const double eh, const double ih){
    const double term = (eh/2.25);
    const double term2 = (eh/2.84) * (eh/2.84);
    const double term3 = (eh/2.02) * (eh/2.02);
    const double Pe = (1.0 + pow(term, 1.2) +  term2*term2*term2) / (1.0 - term3*term3);
    const double DeltaL = 1.0/(Pe + Pe/fabs(Pe)*(0.070*ih + 0.085*ih*ih*ih*ih - 0.080*eh*ih*ih));
    const double GammaL = (-2.5 - 1.7*temp_ind + 0.1*sd_ind)/adi_ind;
    return DeltaL*GammaL;
}

// Function governing the saturation of the corotation torque, PBK11 Eq. 23
const double F_sat(const double p){
    return 1.0/(1.0 + (p*p/(1.3*1.3)));
}

// Function governing the cutoff of corotation torque from diffusion, PBK11 Eq. 30
const double G_cutoff(const double p){
    if (p < sqrt(8.0/(45.0*M_PI))){
        return 16.0/25.0*pow(45*M_PI/8.0, 0.75)*pow(p, 1.5);
    }
    else{
        return 1.0 - 9.0/25.0*pow(8.0/(45.0*M_PI), 4.0/3.0)*pow(p, -8.0/3.0);
    }
}

// Function governing the cutoff of corotation torque from diffusion, PBK11 Eq. 31
const double K_cutoff(const double p){
    if (p < sqrt(28.0/(45.0*M_PI))){
        return 16.0/25.0*pow(45.0*M_PI/28.0, 0.75)*pow(p, 1.5);
    }
    else{
        return 1.0 - 9.0/25.0*pow(28.0/(45.0*M_PI), 4.0/3.0)*pow(p, -8.0/3.0);
    }
}

const double rebx_calculate_corotation_torque(const double sd_ind, const double temp_ind, const double adi_ind, const double e, const double i, const double h, const double alpha_visc, const double mass_ratio){
    const double xs = 1.1/pow(adi_ind, 0.25) * sqrt(mass_ratio/h);
    const double p_nu = 2.0/3.0*sqrt(xs*xs*xs/(2.0 * M_PI * alpha_visc * h * h));
    const double p_chi = 1.5*p_nu;
    const double F_p_nu = F_sat(p_nu);
    const double F_p_chi = F_sat(p_chi);
    const double G_p_nu = G_cutoff(p_nu);
    const double G_p_chi = G_cutoff(p_chi);
    const double K_p_nu = K_cutoff(p_nu);
    const double K_p_chi = K_cutoff(p_chi);

    const double Gammahsbaro = 1.1*(1.5 - sd_ind)/adi_ind;
    const double Gammalinbaro = 0.7*(1.5 - sd_ind)/adi_ind;
    const double Gammabaro = Gammahsbaro * F_p_nu * G_p_nu + Gammalinbaro * (1.0 - K_p_nu);

    const double ent_ind = temp_ind - (adi_ind - 1.0)*sd_ind;
    const double Gammahsent = 7.9*ent_ind/(adi_ind * adi_ind);
    const double Gammalinent = (2.2 - 1.4/adi_ind) * ent_ind / adi_ind;
    const double Gammaent = Gammahsent * F_p_nu * F_p_chi * sqrt(G_p_nu * G_p_chi) + Gammalinent * sqrt((1.0 - K_p_nu) * (1.0 - K_p_chi));
    const double ef = 0.5*h + 0.01;
    const double DeltaC = exp(-e/ef);
    return DeltaC*(Gammabaro + Gammaent);
}

const double rebx_calculate_total_torque(const double sd_ind, const double temp_ind, const double adi_ind, const double e, const double i, const double h, const double alpha_visc, const double mass_ratio){
    const double GammaL = rebx_calculate_lindblad_torque(sd_ind, temp_ind, adi_ind, e/h, i/h);
    const double GammaC = rebx_calculate_corotation_torque(sd_ind, temp_ind, adi_ind, e, i, h, alpha_visc, mass_ratio);
    return GammaL + GammaC;
    // note--these torques are normalized by the reference torque Gamma0
}

// Total migration timescale
const double rebx_calculate_type_I_migration_timescale(const double wave, const double sd_ind, const double temp_ind, const double adi_ind, const double e, const double i, const double h, const double alpha_visc, const double mass_ratio){
    const double Gamma = rebx_calculate_total_torque(sd_ind, temp_ind, adi_ind, e, i, h, alpha_visc, mass_ratio);
    const double t_mig = wave/h/h/Gamma;
    return t_mig;
}

// Eccentricity damping timescale from Cresswell & Nelson 2008. 
const double rebx_calculate_ecc_damping_timescale(const double wave, const double eh, const double ih){
    const double t_e = (wave/0.780) * (1.0 - (0.14*eh*eh) + (0.06*eh*eh*eh) + (0.18*eh*ih*ih));
    return t_e;
}

// Inclination damping timescale from Cresswell & Nelson 2008
const double rebx_calculate_inc_damping_timescale(const double wave, const double eh, const double ih){
    const double t_i = (wave/0.544) * (1.0 - (0.30*ih*ih) + (0.24*ih*ih*ih) + (0.14*eh*eh*ih));
    return t_i;
}

static struct reb_vec3d rebx_calculate_modify_orbits_with_all_type_I_torques(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){
    double beta;
    double h0;
    double sd0;
    double background_sd_ind;
    double inner_edge_pos = 0.0;
    double inner_edge_width = INFINITY;
    double bumppos = 1.0;
    double bumpwidth = 1.0;
    double bumpheight = 0.0;
    double alpha_visc = 0.0;
    double adi_ind = 1.4;

    const double* const inner_edge_pos_ptr = rebx_get_param(sim->extras, force->ap, "ide_position");
    const double* const inner_edge_width_ptr = rebx_get_param(sim->extras, force->ap, "ide_width");
    const double* const sd0_ptr = rebx_get_param(sim->extras, force->ap, "tIm_surface_density_1");
    const double* const background_sd_ind_ptr = rebx_get_param(sim->extras, force->ap, "tIm_surface_density_exponent");
    const double* const h0_ptr = rebx_get_param(sim->extras, force->ap, "tIm_scale_height_1");
    const double* const beta_ptr = rebx_get_param(sim->extras, force->ap, "tIm_flaring_index");
    const double* const bumppos_ptr = rebx_get_param(sim->extras, force->ap, "tIm_bump_position");
    const double* const bumpwidth_ptr = rebx_get_param(sim->extras, force->ap, "tIm_bump_width");
    const double* const bumpheight_ptr = rebx_get_param(sim->extras, force->ap, "tIm_bump_height");
    const double* const alpha_visc_ptr = rebx_get_param(sim->extras, force->ap, "tIm_alpha_visc");
    const double* const mmax_ptr = rebx_get_param(sim->extras, force->ap, "tIm_max_mass");

    int err=0;
    struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *source, &err);
  
    const double a0 = o.a;
    const double e0 = o.e;
    const double inc0 = o.inc;
    const double mp = p->m;  
    const double ms = source->m;

    if (mmax_ptr != NULL && mp > *mmax_ptr){
        struct reb_vec3d a = {0};
        return a;
    }

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x-source->x;
    const double dy = p->y-source->y;
    const double dz = p->z-source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;
    const double r = sqrt(r2);

    if (beta_ptr != NULL) beta = *beta_ptr;
    if (background_sd_ind_ptr != NULL) background_sd_ind = *background_sd_ind_ptr;
    if (sd0_ptr != NULL) sd0 = *sd0_ptr;
    if (h0_ptr != NULL) h0 = *h0_ptr;
    if (inner_edge_pos_ptr != NULL) inner_edge_pos = *inner_edge_pos_ptr;
    if (inner_edge_width_ptr != NULL) inner_edge_width = *inner_edge_width_ptr;
    if (bumppos_ptr != NULL) bumppos = *bumppos_ptr;
    if (bumpwidth_ptr != NULL) bumpwidth = *bumpwidth_ptr;
    if (bumpheight_ptr != NULL) bumpheight = *bumpheight_ptr;
    if (alpha_visc_ptr != NULL) alpha_visc = *alpha_visc_ptr;

    const double h = h0 * pow(r, beta); 
    const double h2 = h*h;

    const double eh = e0/h;
    const double ih = inc0/h;

    const double G = sim->G;
    const double sd = rebx_calculate_disk_surface_density(sd0, r, background_sd_ind, bumpheight, bumppos, bumpwidth, inner_edge_pos, inner_edge_width);
    double sd_ind = rebx_calculate_disk_surface_density_index(r, background_sd_ind, bumpheight, bumppos, bumpwidth, inner_edge_pos, inner_edge_width);
    // avoid unphyiscally large values of sd_ind
    // might want to remove this later
    // if (sd_ind < -10.0) sd_ind = -10.0;
    // if (sd_ind > 10.0) sd_ind = 10.0;
    const double temp_ind = 1.0 - 2.0*beta;
    const double wave = rebx_calculate_wave_timescale(G, sd, r, ms, mp, a0, h2);
    const double invtau_mig = 1.0/rebx_calculate_type_I_migration_timescale(wave, sd_ind, temp_ind, adi_ind, e0, inc0, h, alpha_visc, mp/ms);
    const double tau_e = rebx_calculate_ecc_damping_timescale(wave, eh, ih);
    const double tau_inc = rebx_calculate_inc_damping_timescale(wave, eh, ih);

    struct reb_vec3d a = {0};

    if (invtau_mig != 0.0){
        a.x = dvx*(invtau_mig);
        a.y = dvy*(invtau_mig);
        a.z = dvz*(invtau_mig);
    }

    if (tau_e < INFINITY || tau_inc < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = -2.0*vdotr/r2/tau_e;
        a.x += prefac*dx;
        a.y += prefac*dy;
        a.z += prefac*dz - 2.0*dvz/tau_inc;
    }
    return a;
}

void rebx_modify_orbits_with_all_type_I_torques(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_with_all_type_I_torques, particles, N);
}
