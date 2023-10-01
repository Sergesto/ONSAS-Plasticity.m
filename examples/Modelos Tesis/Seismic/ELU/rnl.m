% Diseño sísmico, análisis no lineal
% Resolución por diferencias finitas
%
% Ing. Sergio A. Merlino Chiozza
%
% PEER NGA STRONG MOTION DATABASE RECORD
% SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, 164 (CDMG STATION 279)                
% ACCELERATION TIME HISTORY IN UNITS OF G
% 4164    0.0100    NPTS, DT
%
clc;
clear;
close all;
g=9.8182;                               % aceleración g m/s^2
fileID = fopen('PUL164.AT2');
data = textscan(fileID, '%f', 'HeaderLines',4);
fclose(fileID);
F = data{1};
M=F.*g;
for j=1:4164
  T(j)=(j-1)*0.01;
end
%
%
% C=2.902758613235145*0.5;              % coeficiente sísmico
C=0.5;                                  % coeficiente sísmico
Tn=0.40;                                % período natural, segundos
wn=2*pi/Tn;                             % frecuencia natural, rad/s
xi=0.05;                                % coeficiente de amortiguamiento
MT=141.8E3;                             % masa total de la estructura, Kg
K=wn^2*MT;                              % rigidez de la estructura N/m
dt=0.001;
np=floor(T(4164)/dt);
a1=xi*wn/dt+1/(dt^2);
a2=2/(dt^2);
a3=xi*wn/dt-1/(dt^2);
%
%
for j=1:np
    Tx(j)=(j-1)*dt;
end
Mx=interp1(T,M,Tx,'linear');
% Condiciones iniciales
mu(1)=0;
mu(2)=-Mx(1)*(wn^2)/(a1*9.8182*C);
f(1)=0;
for j=2:(np-1)
    f(j)=f1(f(j-1),mu(j)-mu(j-1));
    mu(j+1)=1/a1*(-Mx(j)/9.8182*wn^2/C+a2*mu(j)+a3*mu(j-1)-wn^2*f(j));
    f(j+1)=f1(f(j),mu(j+1)-mu(j));
end
fo=C*MT*g;
fy=fo;
yf=fy/K;
y=mu.*yf;
    figure(1)
        plot(Tx,Mx);
          title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL')
          xlabel('Tiempo cada 0,01s)');
          ylabel('Aceleración m/s^2');
          legend('Üg(t)')
     figure(2)
        plot(mu,f,'r');
          title('f(\mu)');
          xlabel('Ductilidad \mu');
          ylabel('f');
          legend('f(\mu)');
          hold on;
     figure(3)
        plot(Tx,mu,'r');
          title('Ductilidad \mu');
          xlabel('t(s)');
          ylabel('Ductilidad \mu');
          legend('\mu');
          hold on;
     figure(4)
        plot(Tx,y,'b');
          title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL');
          xlabel('Tiempo t');
          ylabel('Desplazamientos y');
          legend(['y(t) Coeficiente sísmico C=',num2str(C)]);
%
disp(sprintf('Ductilidad máxima µ=%0.2f',max(abs(mu))));
%
fy=C*MT*g;
yf=fy/K;
y=mu.*yf;
        figure(4)
        plot(Tx,y,'b');
          title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL');
          xlabel('Tiempo t');
          ylabel('Desplazamientos y');
          legend(['y(t) Coeficiente sísmico C=',num2str(C)]);