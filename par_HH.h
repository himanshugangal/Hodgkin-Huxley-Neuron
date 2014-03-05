
double	freq_min, freq_max, dfreq ;
long   Num_I0, Num_freq;
double Iext0_min, Iext0_max, dIext0;
double t_max;




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
	Num_freq = 1; // number of different frequency values (of current) for which simulations will be carried
	printf("Number of different values of forcing frequencies= %d \n \n",Num_freq);

	// Maximum time for which simulation will be run
	t_max = 100;  // in seconds
	printf("Maximum time(in seconds) = %3.10f \n \n",t_max);


