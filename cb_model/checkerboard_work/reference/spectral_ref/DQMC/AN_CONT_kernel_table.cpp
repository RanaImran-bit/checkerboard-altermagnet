#include "AN_CONT/general_parameters.h"


double kernel_conductivity(double tau, double omega)
{
   if(fabs(omega) < 1.0e-10) {
                // Limit as omega->0
                return (1.0/__PI) * (1.0 / (beta / 2.0));
            } else {
                return (omega / __PI) * 
                    (cosh(omega * (tau - beta / 2.0)) / sinh(omega * beta / 2.0));
            }  
}
double kernel_SP_spectral(double tau, double omega)
{
    if(flag_symmetrization)
			return (1.0/__PI)*(cosh(omega*(tau-beta/2.0))/cosh(omega*beta/2.0));
	else
			return (1.0/(2.0*__PI))*(exp(-omega*tau)/(exp(-omega*beta)+1.0));
}



// OPTIMIZATION: Kernel lookup table
const int N_kernel_omega_points = 50000;  // Fine grid for accurate interpolation

kernel_table::kernel_table()
{
    kernel_table_conductivity =NULL;
    kernel_omega_min=0.0;
    kernel_omega_max=0.0;
    kernel_domega=0.0;
}

void kernel_table::format(double (*kernel_func_in)(double, double), FILE* log_file)
{
    kernel_func_ptr=kernel_func_in;
    // Allocate table
    kernel_table_conductivity = new double*[Nt_intervals_half];
    for(int t = 0; t < Nt_intervals_half; t++) {
        kernel_table_conductivity[t] = new double[N_kernel_omega_points];
    }
    kernel_omega_min = flag_symmetrization ? 0.0 : -fabs(omega_int_limit_SP) ;
    kernel_omega_max = fabs(omega_int_limit_SP);   
     kernel_domega = (kernel_omega_max - kernel_omega_min) / (double)(N_kernel_omega_points - 1);

     // Pre-compute all kernel values
    for(int t = 0; t < Nt_intervals_half; t++) {
        double tau = euclidean_time(t);
        for(int w = 0; w < N_kernel_omega_points; w++) {
            double omega = kernel_omega_min + w * kernel_domega;
            
            // Compute kernel K(tau, omega) for conductivity
             kernel_table_conductivity[t][w]=kernel_func_ptr(euclidean_time(t), omega);
        }
    }
    
   fprintf(log_file, "Kernel table initialized: %d tau points x %d omega points\n", 
           Nt_intervals_half, N_kernel_omega_points);
    fprintf(log_file, "Omega range: [%.3f, %.3f], spacing: %.6f\n", 
           kernel_omega_min, kernel_omega_max, kernel_domega);
    fflush(log_file);

}
kernel_table::~kernel_table()
{
    if(kernel_table_conductivity != NULL) {
        for(int t = 0; t < Nt_intervals_half; t++) {
            delete[] kernel_table_conductivity[t];
        }
        delete[] kernel_table_conductivity;
        kernel_table_conductivity = NULL;
    }
}

 double kernel_table::kernel_lookup(int tau_index, double omega) {
    if(kernel_table_conductivity==NULL)
        return 0.0;
    // Bounds checking
    if(omega < kernel_omega_min || omega > kernel_omega_max) {
        // Fall back to direct computation for out-of-bounds
        return kernel_func_ptr(euclidean_time(tau_index), omega);
    }
    
    // Linear interpolation in table
    double omega_index = (omega - kernel_omega_min) / kernel_domega;
    int idx = (int)omega_index;
    double frac = omega_index - (double) idx;
    
    if(idx >= N_kernel_omega_points - 1 || idx<0) {
         return kernel_func_ptr(euclidean_time(tau_index), omega);
    }
    
    // Linear interpolation between idx and idx+1
    return kernel_table_conductivity[tau_index][idx] * (1.0 - frac) + 
           kernel_table_conductivity[tau_index][idx + 1] * frac;
}



