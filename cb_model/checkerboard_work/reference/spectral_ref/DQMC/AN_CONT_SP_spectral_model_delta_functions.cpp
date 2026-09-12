#include "AN_CONT/single_particle_spectral.h"


	SP_model_delta_functions::SP_model_delta_functions()
	{
		conductivity_parameters=new double[number_parameters];
		conductivity_initial_parameters=new double[number_parameters];
		indexes_varying_parameters=new int[number_varying_parameters];
		min_values=new double[number_varying_parameters];
		max_values=new double[number_varying_parameters];

		GK_kernel_ptr=&global_GK_kernel;
		
		for(int count_parameter=0; count_parameter<number_parameters; count_parameter++)
		{
		    indexes_varying_parameters[count_parameter]=count_parameter;
		    //even indexes - positions of delta functions, odd - coefficients in fron of them
		    if(count_parameter%2==0)
			{
				if(flag_symmetrization)//if correlator is symetrized - only positive frequencies are considered, otherwise yll frequencies are taken into account
					conductivity_initial_parameters[count_parameter]= omega_int_limit_SP/((double)number_parameters+1.0)*((double)count_parameter+1.0);
				else
					conductivity_initial_parameters[count_parameter]= -omega_int_limit_SP + 2.0*omega_int_limit_SP/((double)number_parameters+1.0)*((double)count_parameter+1.0);
			}
			else
            {
                if(!force_sum_rule)
			        conductivity_initial_parameters[count_parameter]=0.01;
                else
                    conductivity_initial_parameters[count_parameter]=constant_sum_rule/(double) N_SP_delta_functions_omega;
            }    
        }

		for(int count_p=0; count_p<number_parameters; count_p++)
		{
			conductivity_parameters[count_p]=conductivity_initial_parameters[count_p];
		}
	}
	SP_model_delta_functions::~SP_model_delta_functions()
	{
		delete[] conductivity_parameters;
		delete[] conductivity_initial_parameters;
		delete[] min_values;
		delete[] max_values;
		delete[] indexes_varying_parameters;
	}
	SP_model_delta_functions::SP_model_delta_functions(const SP_model_delta_functions &A)
	{
		number_parameters=A.number_parameters;
		number_varying_parameters=A.number_varying_parameters;

		GK_kernel_ptr=A.GK_kernel_ptr;

		conductivity_parameters=new double[number_parameters];
		conductivity_initial_parameters=new double[number_parameters];
		indexes_varying_parameters=new int[number_varying_parameters];
		min_values=new double[number_varying_parameters];
		max_values=new double[number_varying_parameters];

		for(int count_p=0; count_p<number_parameters; count_p++)
		{
			conductivity_parameters[count_p]=A.conductivity_parameters[count_p];
			conductivity_initial_parameters[count_p]=A.conductivity_initial_parameters[count_p];
			if(count_p<number_varying_parameters)
			{
				indexes_varying_parameters[count_p]=A.indexes_varying_parameters[count_p];
				min_values[count_p]=A.min_values[count_p];
				max_values[count_p]=A.min_values[count_p];
			}
		}
	}
	SP_model_delta_functions& SP_model_delta_functions::operator = (const SP_model_delta_functions &A)
	{

		delete[] conductivity_parameters;
		delete[] conductivity_initial_parameters;
		delete[] min_values;
		delete[] max_values;
		delete[] indexes_varying_parameters;

		number_parameters=A.number_parameters;
		number_varying_parameters=A.number_varying_parameters;


		GK_kernel_ptr=A.GK_kernel_ptr;

		conductivity_parameters=new double[number_parameters];
		conductivity_initial_parameters=new double[number_parameters];
		indexes_varying_parameters=new int[number_varying_parameters];
		min_values=new double[number_varying_parameters];
		max_values=new double[number_varying_parameters];

		for(int count_p=0; count_p<number_parameters; count_p++)
		{
			conductivity_parameters[count_p]=A.conductivity_parameters[count_p];
			conductivity_initial_parameters[count_p]=A.conductivity_initial_parameters[count_p];
			if(count_p<number_varying_parameters)
			{
				indexes_varying_parameters[count_p]=A.indexes_varying_parameters[count_p];
				min_values[count_p]=A.min_values[count_p];
				max_values[count_p]=A.min_values[count_p];
			}
		}

		return *this;

	}

//to change for a new model
	double SP_model_delta_functions::width_proposal(int count_var_par)
	{
	    if(indexes_varying_parameters[count_var_par]%2==0)
	    {
		return omega_int_limit_SP/4.0;
	    }
	    else
	    {
		return 0.1;
	    }
    
	}


//to change for a new model
// OPTIMIZED VERSION: Uses kernel table lookup
	void SP_model_delta_functions::corr_from_spectral_func_model(double* current_correlator)
	{
		for(int count_euclidean_time=0; count_euclidean_time<Nt_intervals_half; count_euclidean_time++)
		{
			current_correlator[count_euclidean_time]=0.0;
			for(int count_parameter=0; count_parameter<number_parameters/2; count_parameter++)
			{
				double current_omega=this->conductivity_parameters[2*count_parameter];
				double current_C=this->conductivity_parameters[2*count_parameter+1];

				// FAST: Table lookup instead of expensive cosh computation
				current_correlator[count_euclidean_time] += current_C * 
				   GK_kernel_ptr->kernel_lookup(count_euclidean_time, current_omega);
			}
		}
	}


	//update of current correlator taking into account only change of parameter index_updated_parameter
	void SP_model_delta_functions::update_current_correlator(double* current_correlator)
	{
		for(int count_euclidean_time=0; count_euclidean_time<Nt_intervals_half; count_euclidean_time++)
		{
			if(index_updated_parameter%2==0)
			{
				double current_omega=this->conductivity_parameters[index_updated_parameter];
				double old_omega=stored_parameter_value;

				double current_C=this->conductivity_parameters[index_updated_parameter+1];

				current_correlator[count_euclidean_time] += current_C * (GK_kernel_ptr->kernel_lookup(count_euclidean_time, current_omega) - GK_kernel_ptr->kernel_lookup(count_euclidean_time, old_omega));
			}
			else
			{
				double current_omega=this->conductivity_parameters[index_updated_parameter-1];

				double current_C=this->conductivity_parameters[index_updated_parameter];
				double old_C=stored_parameter_value;

				current_correlator[count_euclidean_time] += (current_C-old_C) * (GK_kernel_ptr->kernel_lookup(count_euclidean_time, current_omega));
			}
		}
	}


void SP_model_delta_functions::update_variable(int index_for_update, unsigned int thread_rnd_seed, double width_proposal)
{
    double lambda1=((double)rand_r(&(thread_rnd_seed)))/((double)RAND_MAX);
    double lambda2=((double)rand_r(&(thread_rnd_seed)))/((double)RAND_MAX);
	if(lambda1<1.0e-10)
		lambda1=1.0e-10;
    double old_value=this->conductivity_parameters[this->indexes_varying_parameters[index_for_update]];
    this->conductivity_parameters[this->indexes_varying_parameters[index_for_update]]+=width_proposal*sqrt(-2.0*log(lambda1))*cos(2.0*__PI*lambda2);
	//printf("%.15le\n", this->conductivity_parameters[this->indexes_varying_parameters[index_for_update]] - old_value);
    double new_value= this->conductivity_parameters[this->indexes_varying_parameters[index_for_update]];

	index_updated_parameter=this->indexes_varying_parameters[index_for_update];
	stored_parameter_value=old_value;

    if(force_sum_rule && this->indexes_varying_parameters[index_for_update]%2==1) 
    {//correlection of all other coefficients in front of delta functions
        double rescaling_alpha=(constant_sum_rule-new_value)/(constant_sum_rule-old_value);
        for(int count_parameter=1; count_parameter<number_parameters; count_parameter+=2)
		{
            if(count_parameter==this->indexes_varying_parameters[index_for_update])
                continue;
            this->conductivity_parameters[count_parameter]*=rescaling_alpha;
        }
    }
}



void SP_model_delta_functions::return_back_variable(int index_for_update)
{
	int parameter_index=0;
	if(index_for_update<0)
		parameter_index=index_updated_parameter;
	else
		parameter_index=this->indexes_varying_parameters[index_for_update];

    double old_value=this->conductivity_parameters[parameter_index];
    this->conductivity_parameters[parameter_index]=stored_parameter_value;
    double new_value= stored_parameter_value;

    if(force_sum_rule && parameter_index%2==1) 
    {//correlection of all other coefficients in front of delta functions
        double rescaling_alpha=(constant_sum_rule-new_value)/(constant_sum_rule-old_value);
        for(int count_parameter=1; count_parameter<number_parameters; count_parameter+=2)
		{
            if(count_parameter==parameter_index)
                continue;
            this->conductivity_parameters[count_parameter]*=rescaling_alpha;
        }
    }
}


double SP_model_delta_functions::control_sum_rule()
{
    double res=0.0;
    for(int count_parameter=1; count_parameter<number_parameters; count_parameter+=2)
	{
        res=res+this->conductivity_parameters[count_parameter];
    }
    return res;
}

//to change for a new model
	bool SP_model_delta_functions::parameter_limits()
	{
		for(int count_parameter=0; count_parameter<number_parameters/2;count_parameter++)
		{
			double current_omega=this->conductivity_parameters[2*count_parameter];
			double current_C=this->conductivity_parameters[2*count_parameter+1];
			if(flag_symmetrization)
			{
				if (current_omega<-1.0e-10)
			   	 return false;
			}
			if (current_C<-1.0e-10)
			    return false;
		}
		return true;
	}

//to change for a new model
	bool SP_model_delta_functions::parameter_limits_local(int count_var_par)
	{
		int count_parameter=this->indexes_varying_parameters[count_var_par]; 
		{
			if(count_parameter%2==0)
			{
			    double current_omega=this->conductivity_parameters[count_parameter];
				if(flag_symmetrization)
				{
					if (current_omega<-1.0e-10)
			   		 return false;
				}

			}
			else
			{
			    double current_C=this->conductivity_parameters[count_parameter];
			    if (current_C<-1.0e-10)
				return false;
			}
		}
		return true;
	}

//to change for a new model
	double SP_model_delta_functions::value(double omega_min, double omega_max)
	{
		double omega=0.5*(omega_min+omega_max);
		
		if(flag_symmetrization)
		{
			if(omega<0.0)
		    	return 0;
		}
		
		double result=0.0;
		
		for(int count_parameter=0; count_parameter<number_parameters/2;count_parameter++)
		{
		    double current_omega=this->conductivity_parameters[2*count_parameter];
		    double current_C=this->conductivity_parameters[2*count_parameter+1];

		    if(current_omega>omega_min && current_omega<omega_max)
		    {
			result=result+current_C;
		    }
		}
		return result/(omega_max-omega_min);
	}

	void SP_model_delta_functions::print(FILE* current_log_file, double* current_correlator)
	{
		fprintf(current_log_file, "\n\nConductivity model output start\n");

		for(int count_par=0; count_par<number_parameters; count_par++)
		{
			fprintf(current_log_file, "%.15le\n", conductivity_parameters[count_par]);
		}
		if(current_correlator!=NULL)
		{
//			fprintf(current_log_file, "action=%.15le\n", deviation_action<SP_model_delta_functions>(current_correlator, this));
			fprintf(current_log_file, "current reconstructed from positive average_sigma (and initial one):\n");
			double* current_correlator_reconstructed=new double[Nt_intervals_half];
			this->corr_from_spectral_func_model(current_correlator_reconstructed);
    		for(int count_j=0; count_j<Nt_intervals_half; count_j++)
    		{
    			fprintf(current_log_file, "%.15le\t%.15le\n", current_correlator_reconstructed[count_j], current_correlator[count_j]);
    		}
    		delete[] current_correlator_reconstructed;
		}
		fprintf(current_log_file, "Conductivity model output end\n\n");fflush(current_log_file);
	}
	void SP_model_delta_functions::print_string(FILE* current_log_file, double* current_correlator)
	{
		for(int count_par=0; count_par<number_parameters; count_par++)
		{
		    fprintf(current_log_file, "%.15le\t", conductivity_parameters[count_par]);
		}
		if(current_correlator!=NULL)
		{
//			fprintf(current_log_file, "\t%.15le\t", deviation_action<SP_model_delta_functions>(current_correlator, this));
			double* current_correlator_reconstructed=new double[Nt_intervals_half];
			this->corr_from_spectral_func_model(current_correlator_reconstructed);
    			for(int count_j=0; count_j<Nt_intervals_half; count_j++)
    			{
    			    fprintf(current_log_file, "\t%.15le\t%.15le\t", current_correlator_reconstructed[count_j], current_correlator[count_j]);
    			}
    		    delete[] current_correlator_reconstructed;
		}
	}


