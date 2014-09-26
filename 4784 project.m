sm= 100; 
step= 0.1;
t=0:STEP:SM

change= [0];
current=[50];

%constant
gk= 36;
gna= 120;
gl= 0.3;
E_K=-12;
E_Na= 115;
E_L= 10.6
C=1
Vrest=-70;

%Equations
alpha_m= 0.1*((25-Vm)/(exp((25-vm)/10)-1));
beta_m= 4*exp(-vm/18);
alpha_n= 0.01*((10-vm)/(exp((10-vm)/10)-1));
beta_n= 0.125*exp(-vm/80);
alpha_h= 0.07*exp(-vm/20);
beta_h=1/(exp((30-vm)/10)+1);

m0= (alpha_m)/(alpha_m + beta_m);
n0= (alpha_n)/(alpha_n + beta_n);
h0= (alpha_h)/(alpha_h + beta_h);

Ina= m^3*gna*h*(vm - E_Na);
IK= n^4*gk*(vm - E_K);
IL= gl*(vm - E_L);
Iion= I - IK- INa- IL

