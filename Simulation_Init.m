%%  1 DOF Simulation Initialization File
clear;  clc
computer = 1;                           %   Computer reference for John
%%  Initialize Simulation Parameters 
Env = Environment;                      %   Environmental parameters 
Veh = Vehicle;                          %   Vehicle parameters 
Mot = F15_Motor;                        %   Motor parameters
%%  Run Simulation
sim('Simulation')
%%  Save Outputs
Log_Data
%%  Plot Results
Plot_Kinematics(Sim);                   %   Plot altitude, velocity, and acceleration 