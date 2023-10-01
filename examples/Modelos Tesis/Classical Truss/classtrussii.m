%md# Classical Truss

close all; clear; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters (N/cm2)
E = 210e3 ;
Kplas = 529.5 ;
sigma_Y_0 = 123.6 ;
Fu = -1 ;

global Rstress Rstrain Rstrainacum Largo ;

% x and z coordinates of node 2
x2 = 220 ;
z2 = 127 ;

materials.hyperElasModel  = 'isotropicHardening_logstrain' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;
elements(3).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs = 3 ;
boundaryConds(2).imposDispVals = 0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (t<=1)*(Fu)*t 
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

initialConds = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

% MEBI [Material Element Boundary_Conditions Initial_Conditions]

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 2 3 ] ;

analysisSettings.deltaT        =       1 ;
analysisSettings.finalTime     =    1000 ;

analysisSettings.stopTolDeltau =   1e-12 ;
analysisSettings.stopTolForces =   1e-12 ;
analysisSettings.stopTolIts    =   15   ;


analysisSettings.posVariableLoadBC = 2 ;

otherParams = struct();
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

otherParams.problemName       = 'Arc-Length_Logarithmic_Strain' ;
analysisSettings.methodName   = 'arcLength' ;
analysisSettings.finalTime    = 1000 ;
analysisSettings.incremArcLen = [0.3*ones(1,1000)] ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.1)/1000 ;
analysisSettings.posVariableLoadBC = 2 ;

global arcLengthFlag
arcLengthFlag = 2 ;

global dominantDofs
dominantDofs = 6+5 ;

global scalingProjection
scalingProjection = -1 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDisps_logarithmic_strain =  -matUs(6+5,:) ;
loadFactors_logarithmic_strain  =  loadFactorsMat(:,2) ;

deltas = [-matUs(6+1,:)' -matUs(6+5,:)'] ;
L(:,:)=[0 0; Largo(:,:)];
x1=[0; Rstrain(:,1)];
y1=[0; Rstress(:,1)];
% Fuerza aplicada vertical
F=2*(z2-deltas(:,2))./L(:,2).*7.*y1;

%{
for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x1(t)=Rstrain(t,1);
        y1(t)=Rstress(t,1);
        % Fuerza aplicada vertical
        F(t)=2*(z2-deltas(t,2))./Largo(t,2).*7.*y1(t);
end
for t = 1:(floor(analysisSettings.finalTime/analysisSettings.deltaT))
        x2(t)=Rstrain(t,2);
        y2(t)=Rstress(t,2);
end
%}
figure(1);
plot(-controlDisps_logarithmic_strain,loadFactors_logarithmic_strain,'LineWidth',2.5);
title('Plasticidad / Barras --ONSAS-- F_{apl}(\delta)');
ylabel('Fuerza Vertical F_{apl}');
xlabel('Desplazamiento \delta');
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
figure(2);
plot(x1,-deltas(1:end,2),'LineWidth',2.5);
title('Plasticidad / Barras --ONSAS-- \sigma(\epsilon)');
ylabel('Desplazamiento \delta');
xlabel('Deformación \epsilon');
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
figure(3);
plot(x1,y1,'LineWidth',2.5);
title('Plasticidad / Barras --ONSAS-- \sigma(\epsilon)');
ylabel('Tensión \sigma');
xlabel('Deformación \epsilon');
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
figure(4);
plot(x1,F,'LineWidth',2);
title('Plasticidad / Barras --ONSAS-- F_{apl}(\delta)');
ylabel('Fuerza Vertical F_{apl}');
xlabel('Deformación \epsilon');
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
figure(5);
plot(-deltas(1:end,2),y1,'LineWidth',2);
title('Plasticidad / Barras --ONSAS-- F_{apl}(\delta)');
ylabel('Tensión \sigma');
xlabel('Desplazamiento \delta');
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

%
% scalar parameters (N/cm2)
E = 210e3;
Kplas = 529.5;
sigma_Y_0 = 123.6;
Ae=7;
x2=220;
z2=127;
%{
% solución analítica
syms Fu ;
t=0 ;
for dw=0:-1:-126
t=t+1;
eqn=log(sqrt((z2+dw)^2+x2^2)/(sqrt(x2^2+z2^2)))==-(sigma_Y_0)/E + (E+Kplas)/(E*Kplas)*((Fu/2/Ae)/(z2+dw)*sqrt((z2+dw)^2+x2^2)+sigma_Y_0);
f(t)=vpasolve(eqn,Fu);
d(t)=dw;
end
figure(1)
plot(d,f,'r',linewidth=2.5);
%}
%
% Solución analítica tramo elástico 1
dw=0:-0.001:-0.3;
Fu=log(sqrt((z2+dw).^2+x2^2)./(sqrt(x2^2+z2^2)))*E*2*Ae.*(z2+dw)./sqrt((z2+dw).^2+x2^2);
figure(1)
plot(dw,Fu,'r',linewidth=2.5);
% Solución analítica tramo plástico 2
dw=-0.3:-1:-127;
Fu=((log(sqrt((z2+dw).^2+x2^2)./(sqrt(x2^2+z2^2))) + (sigma_Y_0)/E)*(E*Kplas)/(E+Kplas)-sigma_Y_0)*2*Ae.*(z2+dw)./sqrt((z2+dw).^2+x2^2);
figure(1)
plot(dw,Fu,'r',linewidth=2.5);

% another material
%
% scalar parameters (N/cm2)
E = 210e3 ;
Kplas = 1093 ;
sigma_Y_0 = 330.3 ;
Fu = -1 ;

% x and z coordinates of node 2
x2 = 220 ;
z2 = 127 ;

materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;
elements(3).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs = 3 ;
boundaryConds(2).imposDispVals = 0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (t<=1)*(Fu)*t 
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

initialConds = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

% MEBI [Material Element Boundary_Conditions Initial_Conditions]

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 2 3 ] ;

analysisSettings.deltaT        =       1 ;
analysisSettings.finalTime     =    1000 ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;


analysisSettings.posVariableLoadBC = 2 ;

otherParams = struct();
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

otherParams.problemName       = 'Arc-Length_Logarithmic_Strain' ;
analysisSettings.methodName   = 'arcLength' ;
analysisSettings.finalTime    = 1000 ;
analysisSettings.incremArcLen = [0.3*ones(1,1000)] ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.1)/1000 ;
analysisSettings.posVariableLoadBC = 2 ;

global arcLengthFlag
arcLengthFlag = 2 ;

global dominantDofs
dominantDofs = 6+5 ;

global scalingProjection
scalingProjection = -1 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDisps_logarithmic_strain =  -matUs(6+5,:) ;
loadFactors_logarithmic_strain  =  loadFactorsMat(:,2) ;

figure(1);
plot(-controlDisps_logarithmic_strain,loadFactors_logarithmic_strain,'LineWidth',1.5);
title('Plasticidad / Barras --ONSAS-- F_{apl}(\delta)');
ylabel('Fuerza Vertical F_{apl}');
xlabel('Desplazamiento \delta');
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%
% another material
%
% scalar parameters (N/cm2)
E = 206e3 ;
Kplas = 1545.1 ;
sigma_Y_0 = 258.3 ;
Fu = -1 ;

% x and z coordinates of node 2
x2 = 220 ;
z2 = 127 ;

materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;
elements(3).elemCrossSecParams = { 'circle' , sqrt(7*4/pi)} ;

boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs = 3 ;
boundaryConds(2).imposDispVals = 0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (t<=1)*(Fu)*t 
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ]  ;

initialConds = struct() ;

mesh.nodesCoords = [   0  0   0 ; ...
                      x2  0  z2 ; ...
                    2*x2  0   0 ] ;

% MEBI [Material Element Boundary_Conditions Initial_Conditions]

mesh.conecCell = cell(5,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 3   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 2   ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 2 3 ] ;

analysisSettings.deltaT        =       1 ;
analysisSettings.finalTime     =    1000 ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;


analysisSettings.posVariableLoadBC = 2 ;

otherParams = struct();
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

otherParams.problemName       = 'Arc-Length_Logarithmic_Strain' ;
analysisSettings.methodName   = 'arcLength' ;
analysisSettings.finalTime    = 1000 ;
analysisSettings.incremArcLen = [0.3*ones(1,1000)] ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.1)/1000 ;
analysisSettings.posVariableLoadBC = 2 ;

global arcLengthFlag
arcLengthFlag = 2 ;

global dominantDofs
dominantDofs = 6+5 ;

global scalingProjection
scalingProjection = -1 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDisps_logarithmic_strain =  -matUs(6+5,:) ;
loadFactors_logarithmic_strain  =  loadFactorsMat(:,2) ;

figure(1);
plot(-controlDisps_logarithmic_strain,loadFactors_logarithmic_strain,'LineWidth',1.5);
title('Plasticidad / Barras --ONSAS-- F_{apl}(\delta)');
ylabel('Fuerza Vertical F_{apl}');
xlabel('Desplazamiento \delta');
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';