clear all;
close all;
clc
% datos
g=9.8182;       % aceleración g m/s^2
m = 141.8e3;    % kg
Tn=0.40;        % período natural s
w=2*pi/Tn;      % frecuencia natural rad/s
k=(w^2)*m;      % rigidez N/m
xi=0.05;        % factor de amortiguamiento
%
% Sismo
%
fileID = fopen('PUL164.AT2');
data = textscan(fileID, '%f', 'HeaderLines',4);
fclose(fileID);
F = data{1};
M=F.*g;
for j=1:4164
  T(j)=j*0.01;
end
  figure
    plot(T,M,'r');
      title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL')
      xlabel('Tiempo (0,01s)');
      ylabel('Aceleración m/s^2');
      legend('Üg(t)')
pause(3);
% transformada rápida de Fourier
% construyo vector de frecuencia
N = length(T);
tfin=T(N);
dw = 2*pi/tfin;
wvec = [0:dw:dw*(N-1)]';
Ag=fft(-M);
%
% Respuesta de desplazamientos (iFFT - método del dominio en frecuencia)
H=1./(w^2-wvec.^2+2*xi.*wvec.*w.*i);
V=H.*Ag;
vv=ifft(V,'symmetric');
figure
    plot(T,vv,'b');
      title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL')
      xlabel('Tiempo (cada 0,01s)');
      ylabel('Desplazamiento  (m)');
    legend('Desplazamientos v(t)');
%
% en régimen elástico C=Sa/g
%
Sd=max(abs(vv))     % Sd desplazamiento espectral
Sa=(w^2)*Sd         % Sa aceleración espectral
C=Sa/g              % Coeficiente sísmico
