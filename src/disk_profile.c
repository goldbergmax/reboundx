#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double rebx_calculate_disk_surface_density(const double sd0, const double r, const double s, const double bump_amp, const double bump_pos, const double bump_width, const double inner_edge, const double inner_edge_width){
    double sd;
    // overall density profile
    sd = sd0*pow(r, -s);
    if (bump_amp != 0.0){
        // Gaussian profile pressure bump
        sd *= (1.0 + bump_amp * exp(-0.5 * pow((r - bump_pos)/bump_width, 2)));
    }
    if (inner_edge != 0.0){
        // inner edge
        sd *= tanh((r - inner_edge)/inner_edge_width);
    }
    return sd;
}

const double rebx_calculate_disk_surface_density_index(const double r, const double s, const double bump_amp, const double bump_pos, const double bump_width, const double inner_edge, const double inner_edge_width){
    double sd_ind;
    // overall density profile
    sd_ind = s;
    if (bump_amp != 0.0){
        // pressure bump
        sd_ind += bump_amp * r * (r - bump_pos)/(bump_width * bump_width * (bump_amp + exp(0.5 * pow((r - bump_pos)/bump_width, 2))));
    }
    if (inner_edge != 0.0){
        // inner edge
        sd_ind += -2.0 * r / inner_edge_width / sinh(2.0 * (r - inner_edge)/inner_edge_width);
    }
    return sd_ind;
}
