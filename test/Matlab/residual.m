function res = residual(R,P,k,R0,n)

res = 8*R^2*(1-sqrt(2)*R0/(sqrt(R0^2+R^2)))*(1-(cos(pi/n))^2)+(1-P*R0/k)*R - R0;