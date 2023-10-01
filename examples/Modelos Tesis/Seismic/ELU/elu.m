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
for i=1:4164
  T(i)=(i-1)*0.01;
end
%
%
C1=2.902758613235145;                   % coeficiente sísmico
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
mu(2)=-Mx(1)*(wn^2)/(a1*9.8182*C1);
f(1)=0;
for j=2:(np-1)
    f(j)=f1(f(j-1),mu(j)-mu(j-1));
    mu(j+1)=1/a1*(-Mx(j)/9.8182*wn^2/C1+a2*mu(j)+a3*mu(j-1)-wn^2*f(j));
    f(j+1)=f1(f(j),mu(j+1)-mu(j));
end
fo=C1*MT*g;
fy=fo;
disp(sprintf('Estructura con fy=%0.2fKN (Régimen Elástico)',fy/1000));
yf=fy/K;
disp(sprintf('Desplazamiento de fluencia yf=%0.2fm',yf));
y=mu.*yf;
    figure(1)
        plot(Tx,Mx);
          title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL')
          xlabel('Tiempo cada 0,01s');
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
        plot(Tx,mu,'b');
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
          legend(['y(t) Régimen Lineal C=',num2str(C1)]);
          hold on
%
% Coeficiente sísmico *0.5
%
C2=2.902758613235145*0.5;               % coeficiente sísmico
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
mu(2)=-Mx(1)*(wn^2)/(a1*9.8182*C2);
f(1)=0;
for j=2:(np-1)
    f(j)=f1(f(j-1),mu(j)-mu(j-1));
    mu(j+1)=1/a1*(-Mx(j)/9.8182*wn^2/C2+a2*mu(j)+a3*mu(j-1)-wn^2*f(j));
    f(j+1)=f1(f(j),mu(j+1)-mu(j));
end
        figure(2)
        plot(mu,f,'b');
          title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL');
          xlabel('Ductilidad \mu');
          ylabel('f');
          legend('f(\mu) Lineal','f(\mu) No Lineal');
        figure(3)
        plot(Tx,mu,'r');
          title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL');
          xlabel('t(s)');
          ylabel('Ductilidad \mu');
          legend('\mu régimen Lineal','\mu régimen No Lineal');
fy=C2*MT*g;
disp(sprintf('Estructura con fy=%0.2fKN (Régimen Plástico)',fy/1000));
yf=fy/K;
disp(sprintf('Desplazamiento de fluencia yf=%0.2fm',yf));
y=mu.*yf;
        figure(4)
        plot(Tx,y,'r');
          title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL');
          xlabel('Tiempo t');
          ylabel('Desplazamientos y');
          legend(['y(t) Régimen Lineal C=',num2str(C1)],['y(t) Régimen No Lineal C=',num2str(C2)]);
%
disp(sprintf('Ductilidad µ=%0.2f',max(abs(mu))));