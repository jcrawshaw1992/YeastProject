function res = lengthresidual(L,P,k,L0,n,y)



% res = P*L*y ...
%     - 2*k*sin(pi/n)*(L-L0) ...
%     - k*L/sin(pi/n)*(1-cos(2*pi/n))*(1-sqrt(L0^2+y^2)/(sqrt(L^2+y^2)));


x = L/2;
z = L/4*(1-cos(pi/n))/sin(pi/n);
x0 = L0/2;
z0 = L0/4*(1-cos(pi/n))/sin(pi/n);

Lhat = sqrt(x^2+y^2+z^2);
L0hat = sqrt(x0^2+y^2+z0^2);

res = P*L*y ...
     - 2*k*sin(2.0*pi/n)*(L-L0) ...
     - 4*k*(Lhat - L0hat)*z/Lhat;