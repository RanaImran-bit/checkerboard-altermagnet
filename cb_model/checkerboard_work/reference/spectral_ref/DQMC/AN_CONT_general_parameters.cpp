#include "AN_CONT/general_parameters.h"


bool flag_symmetrization=true;//symmetrization of the spectral function, if symmetrization is made, the kernel should reflect this too
bool flag_ignore_zero=false;//in this case, the value of correlator at zero Euclidean time is simply ignored


//default values - real values are being read from command line
double beta=2.4;
double delta_euclidean_time=0.1;
int Nt_intervals_half=13; //initially 24 timeslices, then symmetrization - 13 points, zero is NOT ignored


  int model_index=10;//single-particle spectral function
  int N_streams_cond=1024;//

 
  double cond_temperature_start=0.001;
  int N_cond_temperature=10;
  double cond_temperature_factor=2.0;

  long int thermalization_length=1000;
  long int analysis_interval=100;
  long int tempering_interval=5;
  long int total_statistics=1000;

  int N_display_omega_points=50;
  double omega_max_display=10.0;

  double factor_multiplication_current=1.0; 

//parameters for tuning of the acceptance rate - adjusted automatically
double minimum_acceptance_rate=0.1;
double maximum_acceptance_rate=0.9;
double mc_proposal_factor=10.0;

int N_streams_binning=20;
bool  flag_detailed_output=false;

unsigned int global_seed=0;


bool force_sum_rule=false;
double constant_sum_rule=0.0;

double relative_error_tolerance=0.3;

bool flag_only_positive_conductivity=true;//enforces conductivity to be positive (is not applicatble to single particle spectral functions)

double omega_int_limit_SP=12.0;
int N_SP_delta_functions_omega=100;

std::vector<double> euclidena_time_array;
double euclidean_time(int count_euclidean_time)
{
  return euclidena_time_array[count_euclidean_time];
 }

kernel_table global_GK_kernel;