% script to calculate steady state ring of chain of springs with Pressure 
% force P radially outwards.

clear 
close all
P=1.0666e4; %3
k=50; %15
R0 = 0.0015;%1.0;
n=10:5:100;
L=0.012

sinpin = sin(pi./n);

R_scaled_length = k*R0/(k-P*R0)*ones(size(n));

res = 3*P*L*y ...
    - 2*k*sin(pi/n)*(L-L0) ...
    - 2*k*L*sqrt(1-cos(pi/n)^2)*(1-sqrt(L0^2+y^2)/(sqrt(L^2+y^2)));
R_scaled_area_1 = (1-sqrt(1-8*P./k.*R0.*R0.*sinpin)).*k./4./P./R0./sinpin;
R_scaled_area_2 = (1+sqrt(1-8*P./k.*R0.*R0.*sinpin)).*k./4./P./R0./sinpin;
R_scaled_area_condition = 8.*P./k.*R0.*R0.*sinpin;

R_unscaled_length = R0./(1-P./2./k./sinpin);
R_unscaled_area_1 = k/2/P*(1-sqrt(1-4*P*R0/k))*ones(size(n));
R_unscaled_area_2 = k/2/P*(1+sqrt(1-4*P*R0/k))*ones(size(n));
R_unscaled_area_condition = 4*P*R0/k*ones(size(n));
 
L_2D = zeros(size(n));

L0 = 2*R0*sin(pi./n);
y = 8*R0./(2*n/3);
for i=1:length(n)
    L_2D(i) = fzero(@(L) lengthresidual(L,P,k,L0(i),n(i),y(i)),L0(i));end
R_2D = L_2D./2./sin(pi./n);
plot(n, R0*ones(size(n)),'k--',n,R_scaled_length,'-ro',n,R_scaled_area_1,'b-+',n,R_scaled_area_2,'b-o',n,R_scaled_area_condition,'b-x',n,R_unscaled_length,'-ko',n,R_unscaled_area_1,'g-+',n,R_unscaled_area_2,'g-o',n,R_unscaled_area_condition,'g-x',n,R_2D,'k:');

legend('R0','scaled length','scaled area min','scaled area max','scaled area condition','unscaled length','unscaled area min','unscaled area max','unscaled area condition','R_2D')

axis([0,110,0,3*R0])

figure

plot(n,L0,'k',n,L_2D,'r')
legend('L0','L')


