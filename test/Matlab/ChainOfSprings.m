% script to calculate steady state length of chain of springs with fore F
% applied at one end.

clear 
close all
F=1.0;
k=30.0;
R0 = 1.0;
n=1:50;

R= F*n/k+R0;

plot(n,R);