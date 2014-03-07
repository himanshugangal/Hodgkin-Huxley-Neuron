%/ THIS IS A CODE FOR 2 Electrically coupled neurons. 
%// This has been tested against an exact analytical solution (called as Test 1).

% ******************************************************************
% Define dV/dt, dm/dt, dh/dt, dn/dt
% ******************************************************************

KV = @(I,V,m,h,n,GNa,GK,GL,VNa,VK,VL,c) ((I - ((GNa*(m^3)*h*(V-VNa)) + (GK*(n^4)*(V-VK))+(GL*(V-VL))))/c);
Km = @(V,m)  ((((0.1*(25-V))/(exp((25-V)/10) - 1))*(1 - m)) - (m*(4*exp(-V/18))));
Kh = @(V,h)  ((0.07*exp(-V/20))*(1-h)) - ((1/(exp((30-V)/10) + 1))*h);
Kn = @(V,n)  ((0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80)));
KIext=@(time,w,Iext0) (-w*Iext0*sin(w*time));

 
 %//*************************************************
%// Fix the value of simulation variables
%//*************************************************
t =  2000;    %// in mill-seconds
dt = 0.05;  %// in milli-seconds
imax = t/dt;
ivmax = 1
Vstep = 0;

%//***********************************************
%// Fix the value of parameters of H-H equations
%//***********************************************
GNa=120 ;
GK = 36;
GL=0.3;
R = 30;
VNa = 115;
VK = -12;
VL=10.5995;
c = 1;
%//I1 = 32;  %//current into neuron 1
Iext0= 3;   %// in micro-amperes
w=0.0638*2*3.14;   %// in rad/milli-second

%//***************** construction of poincare map for forced neuron ************
%**********************************************************************************
tpncr=2*3.14/w; %//poincare section is taken after a time interval tsecn = 1/Time period of forcing.
npncr=tpncr/dt; %// number of time steps after which poincare section is taken

%//***********************************************
%// Initialize all matrices to zero
%//***********************************************
%// Neuron 1 variables 
%V1 = zeros(0, imax);
%m1 = zeros(0, imax);
%n1 = zeros(0, imax);
%h1 = zeros(0, imax);
V1 = zeros(1,imax+1);
m1 = zeros(1,imax+1);
n1 = zeros(1,imax+1);
h1 = zeros(1,imax+1);

%// Neuron 2 variables 
V2 = zeros(0, imax);
m2 = zeros(0, imax);
n2 = zeros(0, imax);
h2 = zeros(0, imax);
I2 = zeros(0, imax);

time=0;


%//***********************************************************
%// Fix the initial value of V
%//***********************************************************
V1(1)=0;
%V1(1)= 30*rand;
V2(1)=0;


%// ****************************************************************************************************
%// Fix the initial values of m, h, n for neuron 1 & neuron 2 (taking initial dm/dt=0, dh/dt=0, dn/dt=0
%// ****************************************************************************************************
alpham1 = (0.1*(25-V1(1)))/(exp((25-V1(1))/10) - 1);
betam1 = 4*exp(-V1(1)/18);
alphah1 = 0.07*exp(-V1(1)/20);
betah1 = 1/(exp((30-V1(1))/10) + 1);
alphan1 = 0.01*(10-V1(1))/(exp((10-V1(1))/10) - 1);
betan1 = 0.125*exp(-V1(1)/80);

m1(1) = (alpham1/(alpham1 + betam1));
h1(1) = (alphah1/(alphah1 + betah1));
n1(1) = alphan1/(alphan1 + betan1);
%m1(1)= rand;
%h1(1)=rand;
%n1(1)=rand;

 
jpncr=0;
for j=1:ivmax   %// initial voltage-steps
    
    V1(1) = V1(1) + Vstep ;
    V2(1) = V2(1) + Vstep ;
    
    for i =1:imax      %// time-step
      time=time+dt;
       
    %// For First Neuron
   %// Iext(i)=Iext0*sin(w*time);
    Iext(i)=Iext0*cos(w*time);
 %//   I2(i) = ((V1(i) - V2(i))/R);
    I2(i) = 0;%// Un-comment for testing against analytical case (Test 1).  // THIS IS PUT EQUAL TO ZERO HERE to ensure no current flow from neuron 1  
             % //to 2 while updating neuron 1 and a DIFFERENT VALUE WHILE UPDATING NEURON 2 to ensure that neuron 2 also has same external current 
               %//as neuron 1.
   %// k11 = dt*Voltage(Iext,V1(i),m1(i),h1(i),n1(i));
    k11 = dt*KV(Iext(i)- I2(i),V1(i),m1(i),h1(i),n1(i),GNa,GK,GL,VNa,VK,VL,c);
    m11 = dt*Km(V1(i),m1(i));
    h11 = dt*Kh(V1(i),h1(i));
    n11 = dt*Kn(V1(i),n1(i));
    Iext11=dt*KIext(time,w,Iext0);
    %//Iext11=0;
    
   
   % // k12 = dt*Voltage(Iext,V1(i),m1(i),h1(i),n1(i));
    k12 = dt*KV(Iext(i)+0.5*Iext11 - I2(i),V1(i)+(0.5*k11),m1(i)+(0.5*m11),h1(i)+(0.5*h11),n1(i)+(0.5*n11),GNa,GK,GL,VNa,VK,VL,c);
    m12 = dt*Km(V1(i)+(0.5*k11),m1(i)+(0.5*m11));
    h12 = dt*Kh(V1(i)+(0.5*k11),h1(i)+(0.5*h11));
    n12 = dt*Kn(V1(i)+(0.5*k11),n1(i)+(0.5*n11));
    Iext12=dt*KIext(time+0.5*dt,w,Iext0);
    %//Iext12=0;
    
    %//k13 = dt*Voltage(Iext,V1(i),m1(i),h1(i),n1(i));
    k13 = dt*KV(Iext(i)+0.5*Iext12 - I2(i),V1(i)+(0.5*k12),m1(i)+(0.5*m12),h1(i)+(0.5*h12),n1(i)+(0.5*n12),GNa,GK,GL,VNa,VK,VL,c);
    m13 = dt*Km(V1(i)+(0.5*k12),m1(i)+(0.5*m12));
    h13 = dt*Kh(V1(i)+(0.5*k12),h1(i)+(0.5*h12));
    n13 = dt*Kn(V1(i)+(0.5*k12),n1(i)+(0.5*n12));
    Iext13=dt*KIext(time+0.5*dt,w,Iext0);
    %//Iext13=0;
    
    %//k14 = dt*Voltage(Iext,V1(i),m1(i),h1(i),n1(i));  
    k14 = dt*KV(Iext(i)+Iext13 - I2(i),V1(i)+k13,m1(i)+m13,h1(i)+h13,n1(i)+n13,GNa,GK,GL,VNa,VK,VL,c);
    m14 = dt*Km(V1(i)+k13,m1(i)+m13);
    h14 = dt*Kh(V1(i)+k13,h1(i)+h13);
    n14 = dt*Kn(V1(i)+k13,n1(i)+n13); 
   
   
    V1(i+1) = V1(i) + ((k11 + 2*k12 + 2*k13 + k14)/6);
    m1(i+1) = m1(i) + ((m11 + 2*m12 + 2*m13 + m14)/6);
    h1(i+1) = h1(i) + ((h11 + 2*h12 + 2*h13 + h14)/6);
    n1(i+1) = n1(i) + ((n11 + 2*n12 + 2*n13 +n14)/6);
   
   
    %//Jion1(i) = GNa*((m1(i))^3)*h1(i)*(V1(i)-VNa) + GK*((n1(i))^4)*(V1(i)-VK);
    %//JNa1(i) = GNa*((m1(i))^3)*h1(i)*(V1(i)-VNa);
    %//JK1(i) = GK*((n1(i))^4)*(V1(i)-VK);
    %//JL(i) = GL*(V1(i)-VL);
    %//Netcurrent1(i)=Iext(i) - I2(i)-JNa1(i)-JK1(i)-JL(i);
   
   
    % If int(i/npncr)*npncr == i then
      % jpncr=jpncr+1;
      % V1pncr(jpncr)=V1(i);
     %end
       
end
%xgrid
plot((1:imax)*dt,V1(1:imax),(1:imax)*dt,Iext(1:imax))
%//plot(1:imax,V1(1:imax))
 
 
end

