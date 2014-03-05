/****************************************************************************
 *                                                                          *
 * File    : main.c                                                         *
 *                                                                          *
 * Purpose : Console mode (command line) program.                           *
 *                                                                          *
 * History : Date      Reason                                               *
 *           00/00/00  Created                                              *
 *                                                                          *
 ****************************************************************************/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "plot.h"


double KV( double I, double V,double m,double h,double n,double GNa,double GK,double GL,double VNa,double VK,double VL, double c );
double Km( double V , double m );
double Kh( double V , double h );
double Kn( double V , double n );
double KIext ( double time , double w , double Iext0 );



int main()
{
		

	FILE *fV_vs_t_id;
	clock_t start, end;
	fV_vs_t_id = fopen("V_vs_t.txt", "w");
		
	FILE *out1;
	out1 = fopen("m.txt", "w");
		
	FILE *out2;
	out2 = fopen("h.txt", "w");

	FILE *out3;
	out3 = fopen("n.txt", "w");
		
	FILE *fI_vs_t_id;
	fI_vs_t_id = fopen("I_vs_t.txt", "w");
		

  	double runTime;
	double dt,T,w,freq;
 //	double t_max;
       // double	freq_min, freq_max, dfreq ;
	double  time;
	double GNa,GK,GL,VNa,VK,VL,c,R;
	double alpham1,betam1,alphah1,betah1,alphan1,betan1 ;
	double k11,m11,n11,h11,Iext11;
	double k12,m12,n12,h12,Iext12;
	double k13,m13,n13,h13,Iext13;
	double k14,m14,n14,h14;
	long  i,j, jf, jI0;
	long  n_max;
	//double Iext0_min, Iext0_max, dIext0;
	double Iext0;
       // long  Num_I0, Num_freq;
	char *forcing_param[30];

	double *V1 ;
	double *m1 ;
	double *h1 ;
	double *n1 ;
	double *Iext;

	
	start = clock();



// ************************ Fix neuron parameters*****************************************************
	GNa=120 ;
	GK = 36;
	GL=0.3;
	R = 30;
	VNa = 115;
	VK = -12;
	VL=10.5995;
	c = 1;
// *********************** Input parameter values (read from file below)*****************************************************
#include</home/himanshu/NeuronSimulations/OneFrequency/Code/WorkingDir/par_HH.h>




// *********************** Start the FREQUENCY loop **************************************************************		
	freq = freq_min;
	for(jf = 0 ; jf < Num_freq  ; jf++) { 

		//freq = freq + dfreq;

		// ************************** Calculate T, dt, n_max and w (all of which depend on freq - loop is over frequency) *********

		// Time period
		T=(1/freq); // in seconds
		printf("Time period (in seconds) = %3.10f \n \n",T);

		// Time step
		dt=(T/1000); // in seconds (time step is always chosen as (1/1000)th of the time period of ext current)
		printf("Time step (in seconds) = %3.10f \n \n",dt);

		// Number of time steps
		n_max = t_max*freq*1000; // n_max = t_max / dt , where dt = T/1000 (1 time period divided into 1000 time intervals), and T = 1/freq
		printf("Number of time steps  = %d \n \n",n_max);

		// Anglular frequency
		 w = freq*2*M_PI; // angular frequency (radians per second)
		printf("Angular frequency (rad/second) = %3.10f \n \n",w);

	
		// **************************** convert parameters to milliseconds (from seconds) *********************************
	
		// Anglular frequency
		w=(w/1000); // angular frequency (radians per millisecond)
		printf("Angular frequency (rad/milli-second) = %3.10f \n \n",w);
	

		// Time step
		dt=dt*1000; // in milliseconds. Factor of 1000 converts dt to milliseconds.
		printf("Time step (in milli-second) = %3.10f \n \n",dt);

	        V1 = (double *)malloc(n_max*sizeof(double));
		m1 = (double *)malloc(n_max*sizeof(double));
		h1 = (double *)malloc(n_max*sizeof(double));
		n1 = (double *)malloc(n_max*sizeof(double));
		Iext = (double *)malloc(n_max*sizeof(double));

	
		// ***************************** Prepare simulation for a fresh value of forcing amplitude ***************************************
		Iext0 = Iext0_min;
		for(jI0 = 0 ; jI0 < Num_I0  ; jI0++) {  // Loop over a range of current (external) amplitudes (Max number of current values = Num_I0)
			
			 printf("jI0= %f",jI0);
		
			 // Iext0 = Iext0 + dIext0;  External current will be updated at the end of the loop


			 //  ************* Initialize time, V[], m1, h1, n1 before starting RK4 *********************************
			 // Bring time back to zero
			 time=0;		 

			 // Bring back voltage and gate variables back to rest state
			 V1[0] = 0; 		 
			 alpham1 = (0.1*(25-V1[0]))/(exp((25 - V1[0])/10) - 1);
			 betam1 = 4*exp(-V1[0]/18);
			 alphah1 = 0.07*exp(-V1[0]/20);
			 betah1 = 1/(exp((30-V1[0])/10) + 1);
			 alphan1 = 0.01*(10-V1[0])/(exp((10-V1[0])/10) - 1);
			 betan1 = 0.125*exp(-V1[0]/80);
			 m1[0] = (alpham1/(alpham1 + betam1));
			 h1[0] = (alphah1/(alphah1 + betah1));
			 n1[0] = (alphan1/(alphan1 + betan1));

		
		         fprintf(fV_vs_t_id,"\n");
		         fprintf(fV_vs_t_id,"\n");
		         fprintf(fV_vs_t_id,"%s %f \n","#Iext0=",Iext0);
			 fprintf(fV_vs_t_id,"%s %f \n","#Freq=",freq);
			 fprintf(fV_vs_t_id,"%d %f %f \n",0,time*0.00,V1[0]);

			 fprintf(fI_vs_t_id,"\n");
		         fprintf(fI_vs_t_id,"\n");
		         fprintf(fI_vs_t_id,"%s %f \n","#Iext0=",Iext0);
			 fprintf(fI_vs_t_id,"%s %f \n","#Freq=",freq);

		         

			 // ***************************** START TIME STEPPING WITH RK4 ***************************************
			 for(i=0; i<n_max ; i++) // Loop over time for a chosen forcing parameter values
			 {
				  //printf("i= %d \n",i);
				 //time = time + dt;
				 time = (i+1)*dt;
				 //Time[i] = time ;

				 Iext[i]=Iext0*cos(w*time);


				 k11 = dt*KV(Iext[i],V1[i],m1[i],h1[i],n1[i],GNa,GK,GL,VNa,VK,VL,c);
	   			 m11 = dt*Km(V1[i],m1[i]);
	    			 h11 = dt*Kh(V1[i],h1[i]);
	    			 n11 = dt*Kn(V1[i],n1[i]);
	    			 Iext11=dt*KIext(time,w,Iext0);

				 k12 = dt*KV(Iext[i]+0.5*Iext11 ,V1[i]+(0.5*k11),m1[i]+(0.5*m11),h1[i]+(0.5*h11),n1[i]+(0.5*n11),GNa,GK,GL,VNa,VK,VL,c);
	    			 m12 = dt*Km(V1[i]+(0.5*k11),m1[i]+(0.5*m11));
	    			 h12 = dt*Kh(V1[i]+(0.5*k11),h1[i]+(0.5*h11));
	   			 n12 = dt*Kn(V1[i]+(0.5*k11),n1[i]+(0.5*n11));
	  			 Iext12=dt*KIext(time+0.5*dt,w,Iext0);

				 k13 = dt*KV(Iext[i]+0.5*Iext12 ,V1[i]+(0.5*k12),m1[i]+(0.5*m12),h1[i]+(0.5*h12),n1[i]+(0.5*n12),GNa,GK,GL,VNa,VK,VL,c);
	    			 m13 = dt*Km(V1[i]+(0.5*k12),m1[i]+(0.5*m12));
	    			 h13 = dt*Kh(V1[i]+(0.5*k12),h1[i]+(0.5*h12));
	   			 n13 = dt*Kn(V1[i]+(0.5*k12),n1[i]+(0.5*n12));
	    			 Iext13=dt*KIext(time+0.5*dt,w,Iext0);

				 k14 = dt*KV(Iext[i]+Iext13 ,V1[i]+k13,m1[i]+m13,h1[i]+h13,n1[i]+n13,GNa,GK,GL,VNa,VK,VL,c);
	    			 m14 = dt*Km(V1[i]+k13,m1[i]+m13);
	   		         h14 = dt*Kh(V1[i]+k13,h1[i]+h13);
	   		         n14 = dt*Kn(V1[i]+k13,n1[i]+n13);


				V1[i+1] = V1[i] + ((k11 + 2*k12 + 2*k13 + k14)/6);
	    			m1[i+1] = m1[i] + ((m11 + 2*m12 + 2*m13 + m14)/6);
	   			h1[i+1] = h1[i] + ((h11 + 2*h12 + 2*h13 + h14)/6);
	   			n1[i+1] = n1[i] + ((n11 + 2*n12 + 2*n13 +n14)/6);
					
				fprintf(fV_vs_t_id,"%d %f %f \n",i+1,time*0.001,V1[i+1]);
				//fprintf(out1," %f %f\n",m1[i+1],time);
				//fprintf(out2," %f %f\n",h1[i+1],time);
				//fprintf(out3," %f %f\n",n1[i+1],time);
				fprintf(fI_vs_t_id," %d %f %f \n",i+1,time*0.001,Iext[i]);
	 					
						 
			} // end of RK4 loop

			Iext0 = Iext0 + dIext0; // change the forcing amplitude value and simulate once again.
		} // end of amplitude loop		
  	
		freq = freq + dfreq; // change the forcing frequency value and simulate once again.
  	} // end of frequency loop

		/*	Plotter p = pl_alloc();
			pl_add(p, n_max, Time, V1, "g*-");
		        pl_plot(p);

		//}
				//fclose(ofp);
				end = clock();
  				runTime = (end - start) / (double) CLOCKS_PER_SEC ;
  				printf ("%f\n", runTime);
                                printf ("%d\n", Num_I0);
  				getchar(); */
	fclose(fV_vs_t_id);
	fclose(out1);
	fclose(out2);
	fclose(out3);
	fclose(fI_vs_t_id); 
	return 0; 
}


	double KV( double I, double V,double m,double h,double n,double GNa,double GK,double GL,double VNa,double VK,double VL, double c )
	{
		double p ;
		p =  ((I - ((GNa*(m*m*m)*h*(V-VNa)) + (GK*(n*n*n*n)*(V-VK))+(GL*(V-VL))))/c);
		return (p) ;
	}

	double Km( double V , double m )
	{
		double q;
		q =  ((((0.1*(25-V))/(exp((25-V)/10) - 1))*(1 - m)) - (m*(4*exp(-V/18))));
		return (q);
	}

	double Kh( double V , double h )
	{
		double r;
		r = ((0.07*exp(-V/20))*(1-h)) - ((1/(exp((30-V)/10) + 1))*h);
		return (r);
	}

	double Kn( double V , double n )
	{
		double s;
		s = (0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80)) ;
		return (s);
	}

	double KIext ( double time , double w , double Iext0 )
	{
		double t;
		t = (-w*Iext0*sin(w*time));
		return (t);

	} 

