Tmax= 100;
dt= 0.01;
t= 0:dt:Tmax;
Time= [0];
Level= [50];
I(1:500)=Level;
I(501:2000) = 0;
I(2001:numel(t))= Level;
gk=36;
gna=120;
gl=0.3;
ek=-12;
ena=115;
el=10.6;
c=1;
V=0;
alpha_n= 0.01 * ((10-V)/(exp(10-V)/10)-1);
beta_n= 0.125*exp(-V/80);
alpha_m= 0.1*((25-V)/(exp((25-V)/10)-1));
beta_m= 4*exp(-V/18);
alpha_h= 0.07*exp(-V/20);
beta_h= 1/(exp((30-V)/10)+1);

n(1) = alpha_n/(alpha_n+beta_n);
m(1) = alpha_m/(alpha_m+beta_m);
h(1) = alpha_h/(alpha_h+beta_h);

for i=1:numel(t)-1
    
    alpha_n(i)= 0.01* ( ( 10-V(i))/(exp((10-V(i))/10)-1));
    beta_n(i) = 0.125*exp(-V(i)/80);
    alpha_m(i) = 0.1*((25-V(i))/(exp((25-V(i))/10)-1));
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = 0.07 * exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
    
    Ina= (m(i)^3)* gna*h(i)*(V(i)-ena);
    Ik= (n(i)^4)*gk*(V(i)-ek);
    Il= gl * (V(i)-el);
    Iion= I(i)-Ik-Ina-Il;
    
    V(i+1)= V(i) + dt*(Iion/c);
    n(i+1) = n (i) + dt*(alpha_n(i)*(1-n(i))-beta_n(i)*n(i));
    m(i+1) = m(i) + dt*(alpha_m(i)*(1-m(i))- beta_m(i)* m(i));
    h(i+1) = h(i) + dt* (alpha_h(i)*(1-h(i)) - beta_h(i)*h(i));
end

V=V-70;

plot(t,V,'Linewidth',3)
hold on
legend({'voltage'});
ylabel('Voltage (mv)')
xlabel('time(ms)');
title('Voltage over Time in simulated')

figure
p1= plot(t,gk*n.^4, 'Linewidth',2);
hold on
p2 = plot(t,gna*(m.^3).*h,'r','LineWidth',2)
legend([p1,p2],'Conductance for Potassium','Conductance for Sodium')
ylabel('Conductance')
xlabel('time(ms)')
title('Conducatnace for K and Na')
