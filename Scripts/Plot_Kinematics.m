function Plot_Kinematics(Sim)
k = 3.2808;
figure()
subplot(3,1,1)                          %   Altitude plot
hold on;    grid on
ylabel('Altitude [ft]')
ax = gca;   ax.FontSize = 12;
plot(Sim.t,Sim.z*k,'b-')
subplot(3,1,2)                          %   Velocity plot
hold on;    grid on
plot(Sim.t,Sim.v*k,'b-')
ylabel('Velocity [ft/s]')
ax = gca;   ax.FontSize = 12;
subplot(3,1,3)                          %   Acceleration plot
hold on;    grid on
plot(Sim.t,Sim.a*k,'b-')
xlabel('Time [s]')
ylabel('Acceleration [ft/$\mathrm{s}^2$]')
ax = gca;   ax.FontSize = 12;
fig = gcf;  fig.Position(4) = 550;  fig.Position(2) = 100;