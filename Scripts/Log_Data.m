%%  Define Output Variable Structure
Sim.t = Env_Out.signals.values(:,1);
Sim.rho = Env_Out.signals.values(:,2);
Sim.v_w = Env_Out.signals.values(:,3);
Sim.z = States_Out.signals.values(:,1);
Sim.v = States_Out.signals.values(:,2);
Sim.a = States_Out.signals.values(:,3);
Sim.F_th = Dist_Out.signals.values(:,1);
Sim.F_g = Dist_Out.signals.values(:,2);
Sim.F_d = Dist_Out.signals.values(:,3);
Sim.F_n = Dist_Out.signals.values(:,4);
Sim.date = datestr(now, 'dd-mmmm_HH-MM');
if computer ==1
    results_path = 'C:\Users\John Jr\Dropbox (UNC Charlotte)\Projects\Rocket\1 DOF Simulation Model\Data Files\Results';
else
    results_path = 'C:\Users\jbore\Dropbox (UNC Charlotte)\Projects\Rocket\1 DOF Simulation Model\Data Files\Results';
end
cond = 1;   test = 0;
while cond == 1
    filename = strcat('1-DOF_Output_',num2str(test),'.mat');
    if exist(filename,'file')
        test = test + 1;
    else
        cond = 0;
    end
end
save(strcat(results_path,'\1-DOF_Output_',num2str(test),'.mat'),'Sim')
