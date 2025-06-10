#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double rebx_calculate_disk_surface_density(const double sd0, const double r, const double s_in, const double s_out, const double bump_pos, const double bump_width){
    const double r_bump = (r - bump_pos)/bump_width;
    double ploc = 0.5*(s_out + s_in) + 0.5*(s_out - s_in) * tanh(r_bump);
    double sd = sd0 * pow(r, -ploc);
    return sd;
}

const double rebx_calculate_disk_surface_density_index(const double sd0, const double r, const double s_in, const double s_out, const double bump_pos, const double bump_width){
    const double r_bump = (r - bump_pos)/bump_width;
    double sd_ind = 0.5*(s_out + s_in) + 0.5*(s_out - s_in) * tanh(r_bump); // NEGATIVE power law slope of the gas surface density profile
    sd_ind += 0.5*(s_out - s_in)*r*log(r)/(cosh(r_bump)*cosh(r_bump)*bump_width);
    return sd_ind;
}
