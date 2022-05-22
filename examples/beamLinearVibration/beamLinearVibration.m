%md# Linear Dynamic Vibration of a Simply Supported Beam
%md
%md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS.m/blob/master/examples/beamLinearVibration/beamLinearVibration.m
%md
%mdIn this tutorial, the dynamic response of a simply supported beam is computed using ONSAS with the linear elastic and corotational formulations. The aim of this example is to validate the numerical implementations using the analytic solution.
%md
%mdThe problem consists in a beam with fixed nodes in both ends. In a selected position of the beam a forced load with time dependency $(F = F_o sin(wt))$ is applied in both perpendicular directions, as it is shown in the figure.
%md
%md```@raw html
%md<img src="../../assets/beamDynamicVibration.svg" alt="structure diagram" width="500"/>
%md```
%md
%mdBefore defining the structs, the workspace is cleaned and the ONSAS directory is added to the path
close all, clear all, addpath( genpath( [ pwd '/../../src'] ) );
%mdExternal forced load parameters, time values, Material and Geometric parameters are defined to find both analytic and numerical solutions.
%md
%mdMaterial scalar parameters
E = 200e9 ; nu = 0.3;  rho = 700;
%mdGeometrical scalar parameters
l = 10 ; ty = .3 ;  tz = .1 ;
Iyy = ty*tz^3/12 ;
Izz = tz*ty^3/12 ;
%md
%md Number of elements
numElements = 10 ;
%md Time and applied forced parameters.
Fo     = 100 ; % N
w      = 2   ; % rad/s
tf     = 8   ; % s
deltat = 0.1 ; % s
%md External Load application node
assert( rem( numElements, 2 ) == 0, 'the number of elements must be even.' )
appNode    = ( numElements ) / 2 + 1       ;
appNodePos = (appNode-1) * l / numElements ;
%md
%md## Analytic solution
%md
%md The dynamic displacement of a forced beam is described by the next differential equation
%md```math
%md EI \frac{\partial^4 w}{\partial x^4}(x,t) + \rho A \frac{\partial^2w}{\partial t^2}(x,t) = f(x,t)
%md```
%md Defining a solution $w(x,t) = W(x)T(t)$ it is possible to find:
%md
%md```math
%md w(x,t) = \frac{2fo}{\rho A l} \sum_{n=1}^{\infty} \frac{1}{w_{n}^2 - w^2} \sin\left(\frac{n \pi a}{l} \right) \sin\left(\frac{n \pi x}{l} \right)\sin(wt)
%md```
%md where $f_0$ is the applied force and $\omega_n$ is the natural frequency.
%md
%md### Numerical computation of the analytic solution
%md
ts = 0:deltat:tf       ; % times vector
xs = 0:l/numElements:l ; % beam mesh
ns = 1:8               ; % number of nodes
%md
%md Natural frecuency mode vibration vector
%md
wnY = ( (ns*pi).^2 ) * sqrt(E*Izz/rho/(ty*tz)/(l^4)) ; % Natural frecuency direction Y
wnZ = ( (ns*pi).^2 ) * sqrt(E*Iyy/rho/(ty*tz)/(l^4)) ; % Natural frecuency direction Z
%md
%md Analytic solution
analyticDisY = 0; analyticDisZ = 0;
for i = 1:length( ns )
  analyticDisY = analyticDisY + (1/(wnY(i)^2 - w^2)) .* sin(i*xs/l*pi) .* sin(i*pi*appNodePos/l) .* sin(w*ts)';
  analyticDisZ = analyticDisZ + (1/(wnZ(i)^2 - w^2)) .* sin(i*xs/l*pi) .* sin(i*pi*appNodePos/l) .* sin(w*ts)';
end
analyticDisY = analyticDisY * (2*Fo/(rho*ty*tz*l) ) ;   analyticDisZ = analyticDisZ * (2*Fo/(rho*ty*tz*l) ) ;
%md
%md## Numerical solution
%md### MEBI parameters
%md
%mdThe modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md### Materials
%md Since the example contains only one rod the fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
%md The first analysis case implements the co-rotational formulation
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu]          ;
materials.density         = rho              ;
%md
%md### Elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a rectangular-cross section with $t_y$ and $t_z$ cross-section dimensions in $y$ and $z$ directions, then the elemTypeGeometry field is:
elements(2).elemCrossSecParams = { 'rectangle' , [ty tz] } ;
elements(2).elemTypeParams     = 1                         ;
%md the consistent mass type is selected to solve the dynamic response of the beam
elements(2).massMatType = 'consistent';
%md
%md### boundaryConds
%md
%md The elements are submitted to two different BC settings. The first BC corresponds to the fixed points
boundaryConds(1).imposDispDofs = [ 1 2 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
%md and the second corresponds to an time dependecy external force
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsTimeFact = @(t) Fo*sin(w*t) ;
boundaryConds(2).loadsBaseVals = [ 0 0 1 0 1 0 ] ;
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
%md### mesh parameters
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%md the following case only differs in the boundary condition and the node number
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  numElements+1 ] ;
%md the following case only differs in the boundary condition and the node number
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  appNode ] ;
%md the beam elements are formed by the first material, the second type of element, and no boundary conditions are applied to any element.
for i=1:numElements,
  mesh.conecCell{ i+3,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   deltat  ;
analysisSettings.finalTime     =   tf   ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;
%md
%md## otherParams
otherParams.problemName = 'coRotationaluniformDynamicBeam';
otherParams.controlDofs = [ appNode 3 ] ;
otherParams.plotsFormat = 'vtk' ;
%md Beam with simple supported nodes in the ends
[coRotMatUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md The second analysis case implements the linear elastic formulation
materials.hyperElasModel  = 'linearElastic' ;
otherParams.problemName = 'linearElasticuniformDynamicBeam';
otherParams.controlDofs = [ appNode 5 ] ;
otherParams.plotsFormat = 'vtk' ;
%md
%md Beam with simple supported nodes in the ends
[linElasMatUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md### Error estimation
dofYendNode = 6*(appNode) - 3;
dofZendNode = 6*(appNode) - 1;
%md time vector
timeVec = linspace( 0, tf, size(coRotMatUs,2) );
%md
%md error estimated for each method in the application node of the external force
diflinearDispUy = linElasMatUs(dofYendNode, :) - analyticDisY(: , appNode);
diflinearDispUz = linElasMatUs(dofZendNode, :) - analyticDisZ(: , appNode);
difcoRotDispUy  = coRotMatUs(dofYendNode, :) - analyticDisY(: , appNode);
difcoRotDispUz  = coRotMatUs(dofZendNode, :) - analyticDisZ(: , appNode);
%md
intAnalyticDisY = (analyticDisY(2:end , appNode) + analyticDisY(1:end-1, appNode))*deltat/2;
intAnalyticDisZ = (analyticDisZ(2:end , appNode) + analyticDisZ(1:end-1 , appNode))*deltat/2;
errlinearDispUy = norm((diflinearDispUy(2:end) + diflinearDispUy(1:end-1))*deltat/2)/norm(intAnalyticDisY);
errcoRotDispUy  = norm((difcoRotDispUy(2:end) + difcoRotDispUy(1:end-1))*deltat/2)/norm(intAnalyticDisY);
errlinearDispUz = norm((diflinearDispUz(2:end) + diflinearDispUz(1:end-1))*deltat/2)/norm(intAnalyticDisZ);
errcoRotDispUz  = norm((difcoRotDispUz(2:end) + difcoRotDispUz(1:end-1))*deltat/2)/norm(intAnalyticDisZ);
%md
%md the numerical resolution is validated for both method and both directions.
verifBoolean =  ( errlinearDispUy <  1e-2 ) ...
             && ( errcoRotDispUy  <  1e-2 ) ...
             && ( errlinearDispUz <  1e-2 ) ...
             && ( errcoRotDispUz  <  1e-2 );
%md
%md Plot parameters:
lw = 2.0 ; lw2 = 1.0 ; ms = 11 ; plotfontsize = 18 ;
%md plot y-axis linear, co-rotational and analytic result 
figure(1), hold on, grid on
plot(timeVec, coRotMatUs(dofYendNode, :),'r-x' , 'linewidth', lw,'markersize',ms )
plot(timeVec, linElasMatUs(dofYendNode, :),'k-o' , 'linewidth', lw,'markersize',ms )
plot(timeVec, analyticDisY(:,appNode),'b' , 'linewidth', lw,'markersize',ms )
legend('coRotational_{disp}','linearElastic_{disp}', 'Analytic_{disp}', 'location', 'eastoutside')
labx = xlabel('time (s)');   laby = ylabel('displacement (m)') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/Uy','-dpng')
%md plot z-axis linear, co-rotational and analytic result 
figure(2), hold on, grid on
plot(timeVec, coRotMatUs(dofZendNode, :),'r-x' , 'linewidth', lw, 'markersize', ms )
plot(timeVec, linElasMatUs(dofZendNode, :),'k-o' , 'linewidth', lw, 'markersize', ms )
plot(timeVec, analyticDisZ(:, appNode), 'b' , 'linewidth', lw, 'markersize', ms )
legend('coRotational_{disp}', 'linearElastic_{disp}', 'Analytic_{disp}', 'location', 'eastoutside')
labx = xlabel('time (s)');   laby = ylabel('displacement (m)') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/Uy','-dpng')
%md