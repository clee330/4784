%ECE 4784
%Changdae Lee
%Project 
%9/29

clear;
clc;


Tmax = 100; % in ms
dt= 0.01; % Step size
t = 0:dt:Tmax; % Time vector

% current
I(1:numel(t)) = 0; %Input current= 0 

% Constants
gk = 36; % mS/c^2
gna = 120; % mS/c^2
gl = 0.3; % mS/c^2
ek = -12; % mV
ena = 115; % mV
el = 10.6; % mV
c = 1.0; % uF/c^2
V = 0; % V
%Initial variabless
%alpha= prob of gate close, beta= prob of gate open.
alpha_n = 0.01*((10-V)/(exp((10-V)/10)-1)); %prob of gate close for k
beta_n = 0.125*exp(-V/80);%prob of gate open for k
alpha_m = 0.1*((25-V)/(exp((25-V)/10)-1));%prob of gate close fo Na
beta_m = 4*exp(-V/18);%prob of gate open for Na
alpha_h = 0.07*exp(-V/20);%prob of gate close for l
beta_h = 1/(exp((30-V)/10)+1);%prob of gate open for l

%Initial Condition
n(1) = alpha_n/(alpha_n + beta_n);%no
m(1) = alpha_m/(alpha_m + beta_m);%mo
h(1) = alpha_h/(alpha_h + beta_h);%ho
for i = 1:numel(t)-1 %Iterate through the number of elements in time vector
% Variables with loop
alpha_m = 0.1*((25-V(i))/(exp((25-V(i))/10)-1));
beta_m = 4*exp(-V(i)/18);
alpha_n = 0.01*((10-V(i))/(exp((10-V(i))/10)-1));
beta_n = 0.125*exp(-V(i)/80);
alpha_h = 0.07*exp(-V(i)/20);
beta_h = 1/(exp((30-V(i))/10)+1);

%Currents with loop
Ina = (m(i)^3)*gna*h(i)*(V(i) - ena); %current for Na
Ik = (n(i)^4)*gk*(V(i)-ek); % Current for K
Il = gl*(V(i)-el); %Current for leakage
Iion = I(i) - Ik - Ina - Il; % Total ion current

%How the variable change with time
dv = Iion/c; %How voltage change with time
dm = alpha_m*(1-m(i)) - beta_m*m(i); % How Na change with time
dn = alpha_n*(1-n(i)) - beta_n*n(i); % How K change with time
dh = alpha_h*(1-h(i)) - beta_h*h(i); % How leakage change with time
% Euler's Method
V(i+1) = V(i) + dt*dv; %Euler's for Voltage
m(i+1) = m(i) + dt*dm; %Euler's for Na
n(i+1) = n(i) + dt*dn; %Euler's for K
h(i+1) = h(i) + dt*dh; %Euler's for I
end
V = V-70; % Because the Vrest is -70mV use that to shift down

subplot(2,1,1); %Use the subplot to give two graph in one window
plot(t,V) %Plot the voltage graph
hold on
legend({'Voltage'}) %naming
ylabel('Voltage (mV)') %label Y axis
xlabel('time (ms)') %label X axis  
axis([0, 100, -100, 100]) %Limit the axis
title('Voltage over Time') %Title

subplot(2,1,2);%Use the subplot to give two graph in one window
p1 = plot(t,gk*n.^4);%Plot the conductance for K graph
hold on
p2 = plot(t,gna*(m.^3).*h,'r');%Plot the conductance for Na graph
legend([p1, p2], 'Conductance for K', 'Conductance for Na')
ylabel('Conductance (mS/c^2)') %label Y axis
xlabel('time (ms)') %label X axis
title('Conductance for K and Na') %Title
