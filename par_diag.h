
float	freq_min, freq_max, dfreq ;
long   Num_I0, Num_freq;
float Iext0_min, Iext0_max, dIext0;
double t_max;
long   n_one_cycle;
float T_start_map_constrn;
float T_start_spike_count;




// *********************** Input parameter values *****************************************************
	// Forcing amplitude
	Iext0_min=1.74;  // least value of external current amplitude
        Iext0_max= 1.74;  // largest value of - do -
        dIext0 = 0;    // Difference between successive values of current amplitude
	Num_I0 = 1; 
	printf("Number of forcing Amplitudes = %d \n \n",Num_I0);

	// Forcing frequency (See HH code)
	freq_min = 50; // (in Hz - per second)
	freq_max=  50; // (in Hz - per second)
	dfreq = 0;      // (in Hz - per second)
	Num_freq =1; // number of different frequency values (of current) for which simulations will be carried
	printf("Number of different values of forcing frequencies= %d \n \n",Num_freq);

	// Maximum time for which simulation will be run
	t_max = 100;  // in seconds
	printf("Maximum time(in seconds) = %3.10f \n \n",t_max);


	// Number of time steps in one cycle of sinusoidal input current (later dt is chosen as dt = T/1000). If a change is made,1000 should be changed both places.
	n_one_cycle=1000;
	printf("Number of time steps in one cycle  = %d \n \n ",n_one_cycle);

	// Adjust dimensions, if required
	long j_mode_lck[100][100];

	// Start Poincare map construction at the following time
	T_start_map_constrn = 10; // time at which poincare map construction starts (in seconds)
	printf("Poincare map constructruction starts at time (in seconds) = %f \n \n ",T_start_map_constrn);

	// Start Poincare map construction & also spike counting at the following time
	T_start_spike_count = 1; // time at which spike counting starts (in seconds)
	printf("Spike counting starts at time (in seconds) = %f \n \n ",T_start_spike_count);
