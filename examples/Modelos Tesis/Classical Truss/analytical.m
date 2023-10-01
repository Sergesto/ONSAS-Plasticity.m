% scalar parameters (N/cm2)
E=206e3;
Kplas=1545.1;
sigma_Y_0=258.3;
Ae=7;
x2=220;
z2=127;
% solución analítica
syms Fu ;
t=0 ;
for dw1=0:-1:-126
t=t+1;
eqn1=log(sqrt((z2+dw1)^2+x2^2)/(sqrt(x2^2+z2^2)))==-(sigma_Y_0)/E + (E+Kplas)/(E*Kplas)*(Fu/2/(z2+dw1)/Ae*sqrt((dw1+z2)^2+(z2)^2)+sigma_Y_0);
f(t)=vpasolve(eqn1,Fu);
d(t)=dw1;
end
plot(d,f,'r',linewidth=2);