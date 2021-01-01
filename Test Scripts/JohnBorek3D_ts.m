%%  1 DOF Simulation Initialization File
clear;  clc
computer = 1;                           %   Computer reference for John
%%  Load Components 
loadComponent('noControl');             %   Controller parameters 
loadComponent('Concord3D');             %   Environment parameters 
loadComponent('SLI2020_SubScale');      %   Vehicle parameters 
%%  Environment Properties 
Env.wind0.setValue([5 0 0]*0.621371,'m/s');     %   Set wind velocity vector at ground level 
Env.windZ.setValue([5 0 0]*0.621371,'m/s');    %   Set wind velocity vector at 5000 ft
%%  Vehicle Properties 
Veh.initPosVecGnd.setValue([0;0;Veh.length.Value-Veh.rCMpad_B.Value(1)],'m');
Veh.initVelVecBdy.setValue([0;0;0],'m/s');
Veh.initEulAng.setValue([0;pi/2;0],'rad');
Veh.initAngVelVec.setValue([0;0;0],'rad/s');
Veh.selectCP.setValue(1,'');
Veh.dragCoef.setValue(0.3,'');
%%  Run Simulation
simWithMonitor('Simulation_Model')
%%  Store and Save Outputs
tsc = signalcontainer(logsout);
%%  Plot Results
% tsc.plotPosition;
% tsc.plotVelocity;
%%  Animate Results 
Veh.animateSim(tsc,.1,'GifTimeStep',.1,'PlotTracer',1==0,'FontSize',12,'Pause',1==1,...
    'ZoomIn',1==1,'SaveGif',1==0,'View',[0 0]);
%%
figure; scl = 3.2808;
subplot(3,2,1); hold on; grid on;
plot(squeeze(tsc.positionVec.Data(1,1,:))*scl,squeeze(tsc.positionVec.Data(3,1,:))*scl,'b-');
xlabel('Position [ft]'); ylabel('Altitude [ft]');
subplot(3,2,2); hold on; grid on;
plot(squeeze(tsc.positionVec.Data(3,1,:))*scl,squeeze(tsc.eulerAngles.Data(2,1,:))*180/pi,'b-');
xlabel('Altitude [ft]'); ylabel('Pitch [deg]');
subplot(3,2,3); hold on; grid on;
plot(squeeze(tsc.positionVec.Data(3,1,:))*scl,squeeze(tsc.FAeroBdy.Data(1,1,:)),'b-');
% plot(squeeze(tsc.positionVec.Data(3,1,:))*scl,squeeze(tsc.FAeroBdy.Data(2,1,:)),'r-');
plot(squeeze(tsc.positionVec.Data(3,1,:))*scl,squeeze(tsc.FAeroBdy.Data(3,1,:)),'g-');
xlabel('Altitude [ft]'); ylabel('$\mathrm{F_{aero}}$ [N]'); legend('X','Z');
subplot(3,2,4); hold on; grid on;
plot(squeeze(tsc.positionVec.Data(3,1,:))*scl,squeeze(tsc.alpha.Data(1,1,:))*180/pi,'b-');
xlabel('Altitude [ft]'); ylabel('alpha [deg]');
subplot(3,2,5); hold on; grid on;
plot(squeeze(tsc.positionVec.Data(3,1,:))*scl,squeeze(tsc.MAeroBdy.Data(2,1,:)),'b-');
xlabel('Altitude [ft]'); ylabel('$\mathrm{M_{aero}}$ [Nm]');
subplot(3,2,6); hold on; grid on;
plot(squeeze(tsc.positionVec.Data(3,1,:))*scl,squeeze(tsc.vAppVec.Data(1,1,:))*scl,'b-');
plot(squeeze(tsc.positionVec.Data(3,1,:))*scl,squeeze(tsc.vAppVec.Data(3,1,:))*scl,'g-');
xlabel('Altitude [ft]'); ylabel('$\mathrm{V_{app}}$ [ft/s]');
set(gcf,'Position',[1.8 41.8 766.4 740.8])
