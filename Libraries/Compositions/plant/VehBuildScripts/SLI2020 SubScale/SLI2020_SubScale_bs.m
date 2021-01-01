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
Veh.motor.position.setValue([20*in2m;0;0],'m');

%%  Vehicle parameters 
Veh.length.setValue(21.375*in2m,'m');
Veh.bodyLength.setValue(18*in2m,'m');
Veh.boatTailLength.setValue(0.875*in2m,'m');
Veh.noseLength.setValue(2.5*in2m,'m');
Veh.diameter.setValue(1.375*in2m,'m');
Veh.baseDiameter.setValue(0.8*in2m,'m');

Veh.rCMpad_B.setValue([11.981*in2m;0;0],'m');
Veh.rCMair_B.setValue([11.5*in2m;0;0],'m');
Veh.rCP_B.setValue([15.1848*in2m;0;0],'m');

Veh.dragCoef.setValue(.3,'');
Veh.sideCoef.setValue(1,'');

Veh.structMass.setValue(.503,'kg');
Ixx = 9.0516198e-05;    Iyy = 4.8957999e-03;    Izz = 4.8958000e-03;
Ixy = 1.1688957e-09;    Ixz = 4.4686242e-09;    Iyz = 0;
Veh.inertia_CM.setValue([Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz],'kg*m^2')
                
Veh.numFins.setValue(4,'');
Veh.finLE.setValue(18.36*in2m,'m');
Veh.finThickness.setValue(0.125*in2m,'m');
Veh.finLength.setValue(1.375*in2m,'m');
Veh.finRootChord.setValue(2*in2m,'m');
Veh.finTipChord.setValue(1*in2m,'m');
Veh.finSweep.setValue(30,'deg');
Veh.finIncidence.setValue(0,'deg');

Veh.noseShape.setValue('Haack','');
Veh.boatLE.setValue(20.5*in2m,'m');
%% save file in its respective directory
saveBuildFile('Veh',mfilename,'variant',["PLANT","SIXDOFDYNAMICS","LIBRARY"]);



