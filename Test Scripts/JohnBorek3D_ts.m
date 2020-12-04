%%  1 DOF Simulation Initialization File
clear;  clc
computer = 1;                           %   Computer reference for John
%%  Load Components 
loadComponent('Concord3D');             %   Environment parameters 
% Veh = Vehicle;                          %   Vehicle parameters 

%%  Environment Properties 
Env.wind0.setValue([5 0 0]*0.621371,'m/s');
Env.windZ.setValue([10 0 0]*0.621371,'m/s');


%%  Run Simulation
simWithMonitor('Simulation_Model')
%%  Store and Save Outputs
tsc = signalcontainer(logsout);
%%  Plot Results
Plot_Kinematics(Sim);                   %   Plot altitude, velocity, and acceleration 