function res = radius_residual(R,P,k,R0,L,n_x,n_y)

sinpin = sin(pi/n_y);
x = L/(n_x -1);

res = 2*P*R*sinpin*x ...
    - 4*k*sinpin*sinpin*(R-R0) ...
    - 2*k*R*(1-cos(2*pi/n_y))*(1-sqrt(((2*R0*sinpin)^2+x^2)/((2*R*sinpin)^2+x^2)));
