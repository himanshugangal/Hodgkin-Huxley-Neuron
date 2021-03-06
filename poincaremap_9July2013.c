#include<stdio.h>
#include<math.h>
#include<stdlib.h>

main(){

	//float Iext0_min, Iext0_max,I0,dIext0;
	float  Iext0;
	//float freq_min, freq_max, dfreq;
	float T, freq, w;
	//float t_max,
	float  dt, time;
	//float T_start_map_constrn ,T_end_map_constrn;
	float T_end_map_constrn;
	//long   Num_I0, Num_freq;
	long   n_max;
	//long   n_one_cycle;
	long   n_start_map_constrn, n_end_map_constrn;
	
	float diff;
	float tjunk;
	long ijunk ;

	long n,i,j,temp, k;
	long kf, kI0;
        long jstart;
	long jmax_map;

	char label1[6], label2[5];
	
	float *V, *Vmap;
	//float V[1000000], Vmap[10000];
	
	FILE *ftimeid, *fV_vs_t_id, *fVmapid, *fjunkid , *fmdlck1, *fmdlck2 , *fVmapdiff;

// *********************** Input parameter values (read from file below) *****************************************************

#include</home/himanshu/NeuronSimulations/OneFrequency/Code/WorkingDir/par_diag.h>




// **************** Open Files ***************************************************************************

	fV_vs_t_id=fopen("V_vs_t.txt","r");
	fVmapid=fopen("VPncrMap.txt","w");
        fjunkid=fopen("Vjunk.txt","w");
	fVmapdiff = fopen("VMapdiff.txt","w");
	fmdlck1 = fopen("ModeLckVsI0.txt","w");
	fmdlck2 = fopen("ModeLckVsFreq.txt","w");



// **************** Initialize variables ******************************************************************
	
	//time=0.;
	//j = 0;


// *********************** START LOOP OVER FORCING FREQUENCY ************************************************

	
	
	for(kf=0 ; kf<Num_freq; kf++){
		
		// Above loop index kf is frequency index. Frequencies are read from a file


		// *********************** START LOOP OVER FORCING AMPLITUDES************************************************
		
		for(kI0 = 0 ; kI0 < Num_I0  ; kI0++) { 
			// Above loop index kI0 is amplitude index. Amplitudes are read from a file

			// frequencies and amplitudes are read first. Frequency is needed to calculate dt and n_max needed later.
			fscanf(fV_vs_t_id,"\n");
			fscanf(fV_vs_t_id,"\n");
			fscanf(fV_vs_t_id,"%s %f \n",label1,&Iext0);
			fscanf(fV_vs_t_id,"%s %f \n",label2,&freq);
			fprintf(fjunkid,"#Iext0 = %f Freq = %f \n \n",Iext0, freq);
			fprintf(fVmapid,"\n \n");
			fprintf(fVmapid,"#Frequency= %f \n",freq);
			fprintf(fVmapid,"#Iext0= %f \n",Iext0);
			printf("******Iext0= %f Frequency= %f *****\n ",Iext0, freq);
			
			//********************************  CALCULATE T, dt, n_max ***********************************************************************
			// ************************dt needed for calculating n_start_map_constrn ********************************************************
			// **************n_max = number of voltage values generated by RK4 (depends on the forcing frequency) ***************************
			

			// Time period
			T=(1/freq); // in seconds
			printf("Time period (in seconds) = %f \n \n",T);

			// Time step
			dt=(T/1000); // in seconds (time step is always chosen as (1/1000)th of the time period of ext current)
			printf("Time step (in seconds) = %f \n \n",dt);

			// Anglular frequency
	 		w = freq*2*M_PI; // angular frequency (radians per second)
			printf("Angular frequency (rad/second) = %f \n \n",w);

			n_max = t_max*freq*1000;  //n_max = t_max/dt; 
			printf("Number of time steps = %d \n \n",n_max);

			// ****************** size allocation to arrays (depends on n_max. n_max depends on frequency)********************* 

			V = (float *)malloc(n_max*sizeof(float));
			Vmap = (float *)malloc(0.01*n_max*sizeof(float));
		
			// ****************************READ VOLTAGE TIME SERIES FROM FILE ( generated by running RK4 in another programme)************
     	   		for(i=0;i<=n_max;i++) 
        		{
			
				fscanf(fV_vs_t_id,"%d %f %f \n",&ijunk,&tjunk,&V[i]); 
				fprintf(fjunkid,"%d %f %f \n",ijunk,tjunk,V[i]); /* this is the print statement in the code that generates the time series */
			} 
		
					
			// ********************************* START CONSTRUCTION OF POINCARE MAP *************************************************************

			// *********************First calculate Poincare map parameters from the ones that were input above **************
			n_start_map_constrn = T_start_map_constrn/dt + 1; 
			n_end_map_constrn = t_max/dt + 1 ; //Once started poincare map construction goes on till the end of the time series

			printf("Poincare map construction starts at timestep  = %d \n \n",n_start_map_constrn);
			printf("Poincare map construction ends at timestep  = %d \n \n",n_end_map_constrn);
			
			i=n_start_map_constrn;
			j=0;
			
			for(i=n_start_map_constrn;i<=n_end_map_constrn;i++)
			{
    		
				time=i*dt;

				if ((i-(i/n_one_cycle)*n_one_cycle) == 0)   
				{  						  // n_one_cycle is the number of time steps in 1 cycle of input current 
		
                   	 	 ++j;
				Vmap[j] = V[i];
				
		    	 	fprintf(fVmapid,"%d %d %6.10f %6.15f \n",i, j,time,Vmap[j]);             
              			} 	       	     
			} 

			jmax_map = j;
     
	       		printf("Constructed poincare map \n");

		// ************ FIND MODE LOCKING RATIOS (time period of response/time period of input = j_mode_lck) *******************************************

       
	        	jstart = 1; // start determining mode locking (from poincare map) when the poincare map index j = jstart
        		j=0;
 			diff=10; /* initialize diff to any large enough value */ 
	
		 	while(fabs(diff) >=0.001 && j<jmax_map)
       			{					  
	    				++j;
					//printf("j= %d \n",j);
	    				diff=Vmap[jstart]-Vmap[jstart+j];	   
	    				fprintf(fVmapdiff,"%d %f %f\n",j,Vmap[jstart+j],fabs(diff));  
			} 
			
			if (j == jmax_map) //if the last value of j in previous loop is the last poincare map index. If that is so, mode locking is not found.
			{
		 	    j_mode_lck[kf][kI0]=0;  // We assign a value zero if no mode locking is found till the end of simulation time
			}
			else 
			{
		 	    j_mode_lck[kf][kI0]=j;
			}
			
			//printf("diff=%3.15f Time period of response=%d\n",fabs(diff),j_mode_lck[kf][kI0]; 

	 		printf("Found mode locking ratio \n");
		} // end of Iext0 loop
	
	
	} // end of frequency loop 	


	// *************************************  SAVE MODE LOCKING RATIOS IN EASY TO PLOT FORMAT ******************************************************

	// Print (on file) mode locking ratio vs Iext0 for each frequency value. This data can be used to plot mode locking ration vs Iext0 for a given freq
	
	freq = freq_min;	
	for(kf=0 ; kf<Num_freq; kf++)
	{
		Iext0=Iext0_min;
		fprintf(fmdlck1,"\n \n");
		fprintf(fmdlck1,"#Frequency= %f \n",freq);
		for(kI0 = 0 ; kI0 < Num_I0  ; kI0++)
		{ 
			fprintf(fmdlck1,"%f %d \n",Iext0,j_mode_lck[kf][kI0]);
			Iext0 = Iext0 + dIext0;
		}
	freq = freq + dfreq;
	}

	// Print (on file) mode locking ratio vs freq for each Iext0 value. This data can be used to plot mode locking ratio vs freq for a given Iext0
	
	Iext0=Iext0_min;
	for(kI0=0 ; kI0<Num_I0; kI0++)
	{
		freq = freq_min	;
		fprintf(fmdlck2,"\n \n");
		fprintf(fmdlck2,"#Iext0= %f \n",Iext0);
		for(kf = 0 ; kf < Num_freq  ; kf++)
		{ 
			fprintf(fmdlck2,"%f %d \n",freq,j_mode_lck[kf][kI0]);
			freq = freq + dfreq;
		}
	Iext0 = Iext0 + dIext0;
	}
        
        	 

//***************Close files **************************************************************************************************************************
	fclose(fV_vs_t_id);
	fclose(fVmapid);
	fclose(fjunkid);
	fclose(fVmapdiff);
	fclose(fmdlck1);
	fclose(fmdlck2);

// **************Code ends **************************************************************************************************************************

}
