%%  Manta Ray Vehicle Build Script for 2 Rotors 
clear; clc
in2m = 0.0254;

PLANT                 = "wind3D";
SIXDOFDYNAMICS        = "sixDoFDynamicsQuat";
LIBRARY               = "SLI2020_SubScale";

%% Create Vehicle Object 
Veh = VEH.vehicle;

%%  Motor properties 
Veh.buildMotor;
Veh.motor.position.setValue([.5;0;0],'m');

%%  Vehicle parameters 
Veh.length.setValue(30*in2m,'m');
Veh.boatTailLength.setValue(0.5*in2m,'m');
Veh.diameter.setValue(1.25*in2m,'m');
Veh.baseDiameter.setValue(1*in2m,'m');

Veh.rCMpad_B.setValue([19*in2m;0;0],'m');
Veh.rCMair_B.setValue([18*in2m;0;0],'m');
Veh.rCP_B.setValue([20*in2m;0;0],'m');

Veh.structMass.setValue(1.5,'kg');
Ixx = 3.3181109e+03;    Iyy = 4.0407857e+03;    Izz = 7.2456248e+03;    
Ixy = 7.6635270e-02;    Ixz = 1.0686861e+02;    Iyz = 0;
Veh.inertia_CM.setValue([Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz],'kg*m^2')
                
Veh.numFins.setValue(4,'');
Veh.finLE.setValue(24*in2m,'m');
Veh.finThickness.setValue(0.125*in2m,'m');
Veh.finLength.setValue(3*in2m,'m');
Veh.finRootChord.setValue(4*in2m,'m');
Veh.finTipChord.setValue(3*in2m,'m');
Veh.finSweep.setValue(30,'deg');
Veh.finIncidence.setValue(0,'deg');

%% save file in its respective directory
saveBuildFile('Veh',mfilename,'variant',["PLANT","SIXDOFDYNAMICS","LIBRARY"]);



