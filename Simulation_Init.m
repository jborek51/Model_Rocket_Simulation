%%  1 DOF Simulation Initialization File
clear;  clc
computer = 1;                           %   Computer reference for John
%%  Initialize Simulation Parameters 
Env = ENV.Environment;                      %   Environmental parameters 
Veh = Vehicle;                          %   Vehicle parameters 
Mot = F15_Motor;                        %   Motor parameters
%%  Run Simulation
simWithMonitor('Simulation_Model')
%%  Store and Save Outputs
tsc = signalcontainer(logsout);
% Log_Data
% Save_Data
%%  Plot Results
Plot_Kinematics(Sim);                   %   Plot altitude, velocity, and acceleration 