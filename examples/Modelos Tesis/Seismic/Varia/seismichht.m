%md# Plasticity | spring-mass system / seismic model              
tic ;
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) ) ;
%
% PEER NGA STRONG MOTION DATABASE RECORD
% SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, 164 (CDMG STATION 279)                
% ACCELERATION TIME HISTORY IN UNITS OF G
% 4164    0.0100    NPTS, DT
%
fileID = fopen('PUL164.AT2') ;
data = textscan(fileID, '%f', 'HeaderLines',4) ;
fclose(fileID);
%
g = 9.8182 ; % aceleration g m/s^2
dt = 0.01 ;
F = data{1} ;
M = F.*g;
for i=1:4164
  T(i)=(i-1)*dt ;
end
figure(1) ;
colororder(["#A2142F"]) ;
plot(T,M,'linewidth', 1.0) ;
title('SAN FERNANDO 02/09/71 14:00, PACOIMA DAM, HORIZONTAL') ;
xlabel('Time step 0,01s') ;
ylabel('Aceleration m/s^2') ;
legend('Üg(t)') ;
%
C1 = 2.902758613235145     ;  % seismic coeficient
Tn = 0.40 ;                   % natural period, seconds
wn = 2*pi/Tn ;                % natural frequency, rad/s
xi = 5/100 ;                  % damping ratio
MT = 141.8e3 ;                % total mass of structure, Kg
K = (wn^2)*MT ;               % stiffness of structure N/m
%
% The following numeric parameters are considered.
% Scalar parameters for spring-mass system / seismic model
%
k    = K ;                  % spring constant
c    = xi*2*MT*wn ;         % damping parameter
m    = MT ;                 % mass of the system
u0   = 0.0 ;                % initial displacement
du0  = 0.0 ;                % initial velocity
%
% Then other parameters are computed:
%
% Numerical solutions
%
% Numerical case 1: truss element model with Newmark method and lumped masses
% 
l   = 100 ;
A   = 1 ;
rho = m * 2 / ( A * l ) ;
E   = k*l/A ;
fo = C1*MT*g/A ;
% elastoplastic parameters (N,m)
Kplas = 0;
sigma_Y_0 = fo ;
% 
% where the material of the truss was selected to set a mass $m$ at the node $2$.
%   
% Materials
% 
materials                    = struct() ;
materials(1).hyperElasModel  = 'isotropicHardening_rotengstrain' ;
materials(1).hyperElasParams = [ E Kplas sigma_Y_0 ] ;
materials(1).density   = 0       ;
materials(2).nodalMass = [m m m] ;
% 
% Elements
% 
% In this case only `'node'` and  `'truss'` elements are considered and the lumped inertial formulation is set for the truss element:
elements             = struct() ;
elements(1).elemType = 'node'                                 ;
elements(2).elemType = 'truss'                                ;
elements(2).elemCrossSecParams = {'circle', [sqrt(4*A/pi) ] } ;
elements(2).massMatType = 'lumped'                            ;
%
% Boundary conditions
%
% The node $1$ is fixed, so the boundary condition set is:
boundaryConds                  = struct() ;
boundaryConds(1).imposDispDofs =  [ 1 3 5 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;
% The node $2$ allows the truss to move in $x$ so the boundary condition set is:
boundaryConds(2).imposDispDofs =  [ 3 5 ] ;
boundaryConds(2).imposDispVals =  [ 0 0 ] ;
% and the external load is added into the same boundary condition using:
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) -MT*M(1+t) ;
boundaryConds(2).loadsBaseVals = [ 1 0 0 0 0 0 ] ;
%
% Initial conditions
% Initial displacement and velocity are set:
aux = zeros(6*2,1) ;  aux(7) = u0 ;
initialConds.U = aux ;
aux(7) = du0 ;
initialConds.Udot = aux ; 
%
% Analysis settings
% alfa HHT method with alfa=0 is equivalent to Newmark
analysisSettings               = struct() ;
% 
analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.alphaHHT      =   0 ;
%
analysisSettings.deltaT        =   1 ;
analysisSettings.finalTime     =   4163 ;
analysisSettings.stopTolDeltau =   1e-10 ;
analysisSettings.stopTolForces =   1e-10 ;
analysisSettings.stopTolIts    =   10 ;
%
% OtherParams
% The nodalDispDamping is added into the model using:
otherParams                  = struct() ;
otherParams.nodalDispDamping = c ;
% The name of the problem is:
%
otherParams.problemName = 'Seismic Model (alfa HHT)' ;
%
% mesh
% Only two nodes are considered so the nodes matrix is:
mesh             = struct() ;
mesh.nodesCoords = [  0  0  0 ; ...
                      l  0  0 ] ;
mesh.conecCell = { } ;
%
% MEBI [Material Element Boundary_Conditions Initial_Conditions]
%
% The first node has no material, the first element of the _elements_ struct, which is `'node'` also the first boundary condition (fixed) and no initial condition is set.
mesh.conecCell{ 1, 1 } = [ 0 1 1    1 ] ;
% The second node has material, the first element of the _elements_ struct, which is `'node'` also the second boundary condition (x disp free) and the first initial condition ($u_0$) is set.
mesh.conecCell{ 2, 1 } = [ 2 1 2    2 ] ;
% Only one element is considered with the first material and the second element setting
mesh.conecCell{ 3, 1 } = [ 1 2 0    1 2 ] ;
%
% Execute ONSAS and save the results:
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
% The numerical displacements of the node $2$ is extracted
vals = matUs(6+1,:) ;
times = linspace(0,analysisSettings.finalTime, size(matUs,2))/100 ;
% Plot verification
%
% The control displacement $u(t)$ is plotted:
figure(2);
hold on, grid on, spanPlot = 8 ; lw = 1.0 ; ms = 11 ; plotfontsize = 11 ;
plot(times(1:spanPlot:end), vals(1:spanPlot:end), 'linewidth', lw)
labx = xlabel('t [s]');   laby = ylabel('u(t) [m]') ;
legend( 'Elastic Model', 'location','northeast')
title('Plasticity | spring-mass system / seismic model') ;
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
hold on;
print('output/seismic model.png', '-dpng' )
%
% /\ PLASTIC MODEL (seismic coeficient multiplied with 0.4) /\
%
C1 = 2.902758613235145*0.4 ;  % seismic coeficient
Tn = 0.40 ;                   % natural period, seconds
wn = 2*pi/Tn ;                % natural frequency, rad/s
xi = 5/100 ;                  % damping ratio
MT = 141.8e3 ;                % total mass of structure, Kg
K = (wn^2)*MT ;               % stiffness of structure N/m
%
% The following numeric parameters are considered.
% Scalar parameters for spring-mass system / seismic model
%
k    = K ;                  % spring constant
c    = xi*2*MT*wn ;         % damping parameter
m    = MT ;                 % mass of the system
u0   = 0.0 ;                % initial displacement
du0  = 0.0 ;                % initial velocity
%
% Then other parameters are computed:
%
% Numerical solutions
%
% Numerical case 1: truss element model with Newmark method and lumped masses
% 
l   = 100 ;
A   = 1 ;
rho = m * 2 / ( A * l ) ;
E   = k*l/A ;
fo = C1*MT*g/A ;
% elastoplastic parameters (N,m)
Kplas = 0;
sigma_Y_0 = fo ;
% 
% where the material of the truss was selected to set a mass $m$ at the node $2$.
%   
% Materials
% 
materials                    = struct() ;
materials(1).hyperElasModel  = 'isotropicHardening_rotengstrain' ;
materials(1).hyperElasParams = [ E Kplas sigma_Y_0 ] ;
materials(1).density   = 0       ;
materials(2).nodalMass = [m m m] ;
% 
% Elements
% 
% In this case only `'node'` and  `'truss'` elements are considered and the lumped inertial formulation is set for the truss element:
elements             = struct() ;
elements(1).elemType = 'node'                                 ;
elements(2).elemType = 'truss'                                ;
elements(2).elemCrossSecParams = {'circle', [sqrt(4*A/pi) ] } ;
elements(2).massMatType = 'lumped'                            ;
%
% Boundary conditions
%
% The node $1$ is fixed, so the boundary condition set is:
boundaryConds                  = struct() ;
boundaryConds(1).imposDispDofs =  [ 1 3 5 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;
% The node $2$ allows the truss to move in $x$ so the boundary condition set is:
boundaryConds(2).imposDispDofs =  [ 3 5 ] ;
boundaryConds(2).imposDispVals =  [ 0 0 ] ;
% and the external load is added into the same boundary condition using:
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) -MT*M(1+t) ;
boundaryConds(2).loadsBaseVals = [ 1 0 0 0 0 0 ] ;
%
% Initial conditions
% Initial displacement and velocity are set:
aux = zeros(6*2,1) ;  aux(7) = u0 ;
initialConds.U = aux ;
aux(7) = du0 ;
initialConds.Udot = aux ; 
%
% Analysis settings
% alfa HHT method with alfa=0 is equivalent to Newmark
analysisSettings               = struct() ;
% 
analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.alphaHHT      =   0 ;
%
analysisSettings.deltaT        =   1 ;
analysisSettings.finalTime     =   4163 ;
analysisSettings.stopTolDeltau =   1e-10 ;
analysisSettings.stopTolForces =   1e-10 ;
analysisSettings.stopTolIts    =   10 ;
%
% OtherParams
% The nodalDispDamping is added into the model using:
otherParams                  = struct() ;
otherParams.nodalDispDamping = c ;
% The name of the problem is:
%
otherParams.problemName = 'Seismic Model (alfa HHT)' ;
%
% mesh
% Only two nodes are considered so the nodes matrix is:
mesh             = struct() ;
mesh.nodesCoords = [  0  0  0 ; ...
                      l  0  0 ] ;
mesh.conecCell = { } ;
%
% MEBI [Material Element Boundary_Conditions Initial_Conditions]
%
% The first node has no material, the first element of the _elements_ struct, which is `'node'` also the first boundary condition (fixed) and no initial condition is set.
mesh.conecCell{ 1, 1 } = [ 0 1 1    1 ] ;
% The second node has material, the first element of the _elements_ struct, which is `'node'` also the second boundary condition (x disp free) and the first initial condition ($u_0$) is set.
mesh.conecCell{ 2, 1 } = [ 2 1 2    2 ] ;
% Only one element is considered with the first material and the second element setting
mesh.conecCell{ 3, 1 } = [ 1 2 0    1 2 ] ;
%
% Execute ONSAS and save the results:
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
% The numerical displacements of the node $2$ is extracted
vals = matUs(6+1,:) ;
times = linspace(0,analysisSettings.finalTime, size(matUs,2))/100 ;
% Plot verification
%
% The control displacement $u(t)$ is plotted:
figure(2) ;
hold on, grid on, spanPlot = 8 ; lw = 1.0 ; ms = 11 ; plotfontsize = 11 ;
plot(times(1:spanPlot:end), vals(1:spanPlot:end), 'linewidth', lw)
labx = xlabel('t [s]');   laby = ylabel('u(t) [m]') ;
legend( 'Elastic Model', 'Plastic Model', 'location','northeast')
title('Plasticity | spring-mass system / seismic model') ;
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/seismic model.png', '-dpng' )
toc