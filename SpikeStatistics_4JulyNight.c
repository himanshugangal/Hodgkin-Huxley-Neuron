#include<stdio.h>
#include<math.h>
#include<stdlib.h>

main(){


float Iext0;
float freq;
float T, dt, w;
float spikegroup;
float AveragingTime,FiringRate;
float tjunk;


long kf, kI0;
long i, j, k;
long m, n, grpnum;
long s, nos;
long ijunk;
long n_max;

long n_start_spike_count, n_end_spike_count;
long Nspikes_tot,Ngrp_tot;
long z;
long ilocmax ;

char label1[6], label2[5];

float *V, *tspike, *tNospike, *ISI;
long *spike_per_grp;


// ***** variables for repetitive sequence
long Start_grpnumL, End_grpnumL, Start_grpnumR, End_grpnumR, diff_grpnumLR;
long imatch;
long grpnumL, grpnumR;
long diff;
long  Last_grp_repeating_seq, Num_grp_repeating_seq ;
long *Spikes_grp_repeating_seq;
long tempL, tempR;
long iexpand, ishift;

// *********



//long spike_per_grp[100000];

FILE *fV_vs_t_id, *fjunkid, *fSpikeTimeid, *fSpike_per_grpid, *fISIid, *fRepeatSeqid;
FILE *fdebug;

// *********************** Input parameter values (read from file below)***************************************************

#include</home/himanshu/NeuronSimulations/OneFrequency/Code/WorkingDir/par_diag.h>



// **************** Open Files ***************************************************************************

	fV_vs_t_id=fopen("V_vs_t.txt","r");
        fjunkid=fopen("Vjunk.txt","w");
	fSpikeTimeid = fopen("SpikeTime.txt","w");
	fSpike_per_grpid = fopen("Spike_per_grp.txt","w");
	fISIid = fopen("ISI.txt","w");
        fRepeatSeqid = fopen("RepeatSeq.txt","w");
	fdebug = fopen("debug.txt","w");



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

			// print labels of frequency and Iext on every file that is created.
			fprintf(fjunkid,"#Iext0 = %f Freq = %f \n \n",Iext0, freq);
			fprintf(fSpikeTimeid,"\n \n");
			fprintf(fSpikeTimeid,"#Iext0 = %f Freq = %f \n",Iext0, freq);
			fprintf(fISIid,"\n \n");
			fprintf(fISIid,"#Iext0 = %f Freq = %f \n",Iext0, freq);
			fprintf(fSpike_per_grpid,"\n \n");
			fprintf(fSpike_per_grpid,"#Iext0 = %f Freq = %f \n",Iext0, freq);
			fprintf(fRepeatSeqid,"\n \n");
			fprintf(fRepeatSeqid,"#Iext0 = %f Freq = %f \n",Iext0, freq);
			


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
			tspike = (float *)malloc((n_max/n_one_cycle)*sizeof(float));
			tNospike = (float *)malloc((n_max/n_one_cycle)*sizeof(float));
			ISI = (float *)malloc((n_max/n_one_cycle)*sizeof(float));
			spike_per_grp = (long *)malloc((n_max/n_one_cycle)*sizeof(long));
			Spikes_grp_repeating_seq = (long *)malloc(n_max*sizeof(long));
		
	
			// *****************************************************************************************************************************
			// ****************************************************************************************************************************
			// *					 READ VOLTAGE TIME SERES FROM A FILE  							*
			// *****************************************************************************************************************************
			// *****************************************************************************************************************************
			
			for(i=0;i<=n_max;i++) 
       			{
			fscanf(fV_vs_t_id,"%d %f %f \n",&ijunk,&tjunk,&V[i]); 
			fprintf(fjunkid,"%d %f %f \n",ijunk,tjunk,V[i]); /* this is the print statement in the code that generates the time series */
			} 


			// *****************************************************************************************************************************
			//*				 CALCULATE - SPIKE TIMINGS & NUMBER OF SPIKES IN EACH GROUP  					*
			// ********************************************************************************************************************************
			

			n_start_spike_count = T_start_spike_count/dt + 1; 
			n_end_spike_count = t_max/dt + 1 ; //Once started poincare map construction goes on till the end of the time series

			printf("Spike counting starts at timestep  = %d \n \n",n_start_spike_count);
			printf("Spike counting ends at timestep  = %d \n \n",n_end_spike_count);

			

			i=n_start_spike_count;
			
			s=0; //n is the spike counter
			nos = 0; // m is the no spike counter


			printf("Number of cycles of input current = %d \n",n_max/n_one_cycle);

			for(grpnum=0;grpnum<=n_max/n_one_cycle;grpnum++) //Initialize arrange spike_count[] = 0. This array stores the number of spike in each group.
			{
			 	spike_per_grp[grpnum]=0;  //grp_spikecount[grpnum] is the number of spikes in the group grpnum
			}

			grpnum = 0;
			//ilocmax = -1;
			for(i=n_start_spike_count;i<=n_end_spike_count-1;i++){
				//printf("i= %d \n",i);
				if(V[i]>V[i+1] && V[i]>V[i-1])
				{
					
					if(V[i]>70)
					{
						s++ ; 
						//printf("s= %d \n",s);                 
						tspike[s]=i*dt ; // tspike[n] is the time at which nth spike (maxima in V) occurs
						//printf("s= %d \n",s);            
						spike_per_grp[grpnum]++ ;   
						ilocmax = +1 ;
						//printf("s= %d \n",s);
						//printf("i= %d grpnum = %d spike_per_grp = %d\n",i, grpnum,spike_per_grp[grpnum]);
						fprintf (fSpikeTimeid,"%d %f %d\n",s,tspike[s],1);
						
					}
					else if(V[i]<70)
					{
						nos++ ;
						tNospike[nos] = i*dt ; // tNospike
						if (ilocmax == +1)
							grpnum++ ; // once a missing spike is detected increment grpnum by one.
						ilocmax = -1;
					
					}

				}
			}
			Ngrp_tot = grpnum; // total number of groups of spikes
			Nspikes_tot = s ; // number of spikes that were generated
			printf("Number of spikes = %d \n",Nspikes_tot);
			printf("Number of spike groups = %d \n",Ngrp_tot);

			//Print - number of spikes in each group
			for(grpnum=1;grpnum<=Ngrp_tot;grpnum++) // grpnum starts from '1' instead of '0'. no spike may appear before a spike. Then grpnum
								// would get incremented from 0 to 1 even befor a spike has appeared.
			{
			   fprintf (fSpike_per_grpid,"%d %d \n",grpnum,spike_per_grp[grpnum]);
			}


			//j=1;

			// *****************************************************************************************************************
			// *		       Calculate & Print -  Interspike interval (ISI) 						   *
			// *****************************************************************************************************************

			for (j=1;j<=Nspikes_tot-1;j++){ // should it be j<=Nspikes_tot-1 or j<Nspikes_tot-1
			  	ISI[j]=tspike[j+1]-tspike[j];
			  	fprintf(fISIid,"%d %f \n",j,ISI[j]);
			}


			// ****************** Histogram of # of spikes in a group ******************************************


			// ****************************************************************************************************************
			// ********************* Average firing rate *****************************************************************
			// ****************************************************************************************************************
			//AveragingTime=(n_end_spike_count - n_start_spike_count)/dt; //.....is this formula correct or does -1 have to be subtracted from top?
			//FiringRate=Nspikes_tot/AveragingTime;
			//fprintf(............)

			
			// ********************************************************************************************************************			
			// ********************************************************************************************************************	
			// *   					IDENTIFY REPETITIVE SPIKE SEQUENCE					       *
			// *********************************************************************************************************************
			// *********************************************************************************************************************





			Start_grpnumL = 1; // The 1st group in all sequences that will be tested for repetivity
					   // (this value will not change hereafter. 1st group in sequence is always group no. 0)
			End_grpnumL = 1;   // the last group in the current sequence that will be tested for repetivity.
					   // if the current sequence is not repetitive, sequence size will be increased by 1.
					   // End_grpnumL will thus be increased by 1.
			iexpand = 0;
			ishift = 0;
			do // this loop expands the size of sequence till a repetitive sequence is found)
			{	
				 iexpand++ ;
				 fprintf(fdebug,"iexpand= %d ishift = %d \n", iexpand, ishift);
				// sequences labeled by L (left in time) will be matched against sequences labeled by R (right in time).
				// Matched 'against' sequence will start at group Start_grpnumR and end at group End_grpnumR.
				// The 1st grp in matched 'against' sequence will start just after the last group of the 'being' matched sequence.

				Start_grpnumR = End_grpnumL+1 ; // the group number on right (in time) where the matching starts right after
								// 'i' is incremented by 1 (that is after R sequence is shifted)
						

				diff_grpnumLR = Start_grpnumR - Start_grpnumL; // difference in group nos where left sequence and right sequence begin.

				imatch = +1;  // if L sequences matches with R sequence , imatch = +1 (imatch also initialized to +1)

				i=1; 
	
        			do   // this loop shifts the R sequence against which the L sequence will be matched.
       				{
					ishift++;
                                        fprintf(fdebug,"ishift= %d iexpand= %d \n", ishift, iexpand);
					// grpnumL is the group label in left(L) sequence. This loop is over groups in a chosen L sequence 
					for(grpnumL=Start_grpnumL; grpnumL<=End_grpnumL; grpnumL++) // loop over all groups in L sequence
					{
	
						grpnumR =  grpnumL+i* diff_grpnumLR ;// (grpnumL will be matched 'against' grpnumR)
									
						//diff = spike_per_grp[grpnumL] - spike_per_grp[grpnumL+i*diff_grpnumLR];
							   // diff in number of spikes in grps being matched in left & right sequences
						tempL = spike_per_grp[grpnumL];
						tempR = spike_per_grp[grpnumL+i*diff_grpnumLR];

					fprintf(fdebug,"grpnumL= %d grpnumR= %d imatch=%d",grpnumL, grpnumR,imatch);						
						fprintf(fdebug,"tempL= %d tempR = %d \n", tempL, tempR);
						
					 	if(tempL != tempR) 											
						//if (fabs(diff) != 0) //(if diff is non-zero sequence L does not repeat)
						{
							imatch = -1; // repetitive sequence not found (imatch = -1)
						}

					}  // loop continues till all groups in sequence L have been compared against the R sequence
					   // (even if imatch == -1 already. not efficient. can be rectified later)

					i++;  // 'i' locates the R sequence. Each R sequence has same length as L sequence.
					      //  But the R sequence can be 1 sequence away, 2 sequences away, or i sequence away from L sequence.

				 }while (imatch == +1 && End_grpnumL+i*diff_grpnumLR <= Ngrp_tot-2);

        			End_grpnumL++ ; // change the number of groups in L sequence by 1 (since the shorter sequence was not repetitive)
				//printf("End_grpnumL= %d \n", End_grpnumL);
				fprintf(fdebug,"Now L sequence will be expanded , imatch= %d \n",imatch);

			}while(imatch == -1 && End_grpnumL <= (Ngrp_tot-2)/2); 
								//if a sequence matches all  other same size sequences imatch=+1 and we
								 //have found a repetitive sequence. Job is done.On the other hand if a 
								//repetitive seq is not found imatch = false and we must go on.
			fprintf(fdebug,"Repetitive sequence found, imatch= %d \n",imatch);
    
			Last_grp_repeating_seq = End_grpnumL -1 ;
			
			Num_grp_repeating_seq = Last_grp_repeating_seq;

			if (imatch == +1)
			{
				for(k = 1;k<= Last_grp_repeating_seq;k++) 
				{
			  		Spikes_grp_repeating_seq[k] = spike_per_grp[Start_grpnumL+k-1];
				 	fprintf(fRepeatSeqid,"%d \n",Spikes_grp_repeating_seq[k]);
				}
			}
			else if (imatch == -1)
			{
				fprintf(fRepeatSeqid,"%d \n", 0);
				printf("No repeating sequence is found");
			}
			



// **********************************************************************************************************************************
		} // amp loop
	} // freq loop
} // main loop










