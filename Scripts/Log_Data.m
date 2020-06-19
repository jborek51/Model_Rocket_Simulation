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
