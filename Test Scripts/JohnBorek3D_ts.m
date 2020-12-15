%%  1 DOF Simulation Initialization File
clear;  clc
computer = 1;                           %   Computer reference for John
%%  Load Components 
loadComponent('noControl');             %   Controller parameters 
loadComponent('Concord3D');             %   Environment parameters 
loadComponent('SLI2020_SubScale');      %   Vehicle parameters 
%%  Environment Properties 
Env.wind0.setValue([5 0 0]*0.621371,'m/s');     %   Set wind velocity vector at ground level 
Env.windZ.setValue([10 0 0]*0.621371,'m/s');    %   Set wind velocity vector at 5000 ft
%%  Vehicle Properties 
Veh.initPosVecGnd.setValue([0;0;Veh.length.Value-Veh.rCMpad_B.Value(1)],'m');
Veh.initVelVecBdy.setValue([0;0;0],'m/s');
Veh.initEulAng.setValue([0;pi/2;0],'rad');
Veh.initAngVelVec.setValue([0;0;0],'rad/s');
Veh.selectCP.setValue(1,'');
%%  Run Simulation
simWithMonitor('Simulation_Model')
%%  Store and Save Outputs
tsc = signalcontainer(logsout);
%%  Plot Results
% tsc.plotPosition;
% tsc.plotVelocity;
%%  Animate Results 
Veh.animateSim(tsc,.1,'GifTimeStep',.1,'PlotTracer',1==0,'FontSize',12,'Pause',1==0,...
    'ZoomIn',1==1,'SaveGif',1==0);
