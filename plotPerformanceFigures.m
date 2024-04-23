%% Create Folder
rootDir = 'Figures\SUAV-I Results';
ARfolder = sprintf('AR %s',num2str(Aircraft.Aero.Wing.AR));
mkdir([rootDir '\' ARfolder]);
resDir = ['Figures\SUAV-I Results\' ARfolder '\' strtok(Mission.dateStart,' ')];
mkdir(resDir);

%% Flight Sim Plotting
time = t./3600;
[~,dayLoc] = min(abs(time - 12));
xp = x(1:dayLoc).*(39.6./12);
yp = y(1:dayLoc).*(39.6./12);
zp = (z(1:dayLoc) - z(1)).*(39.6./12);

xto = xp(takeoffMode(1:dayLoc));
yto = yp(takeoffMode(1:dayLoc));
zto = zp(takeoffMode(1:dayLoc));

zmidClimb = (max(zp(chargingCruiseMode(1:dayLoc))) + min(zp(chargingCruiseMode(1:dayLoc)))) ./ 2;

chLoc = chargingCruiseMode(1:dayLoc) & (zp <= zmidClimb);
crLoc = chargingCruiseMode(1:dayLoc) & (zp > zmidClimb);

xch = xp(chLoc);
ych = yp(chLoc);
zch = zp(chLoc);

xcr = xp(crLoc);
ycr = yp(crLoc);
zcr = zp(crLoc);

xcl = xp(neutralClimbMode(1:dayLoc));
ycl = yp(neutralClimbMode(1:dayLoc));
zcl = zp(neutralClimbMode(1:dayLoc));

% Power and Energy Profile
Cp = (C(1:end-1)./Cpack).*1e2;

figure
subplot(2,1,1)
plot(time,Cp,'linewidth',3)
set(gca,'fontsize',14)
ylabel('Battery Pack Charge (%)')
ylim([0,100.5])
xlim([0,max(time)])
grid on
title(sprintf('SUAV Energy and Power Profile\nLaunch: %s | %d Hour Flight | [%s,%s]',...
    Mission.dateStart,Mission.flightHours,...
    num2str(Mission.location(1)),num2str(Mission.location(2))))

subplot(2,1,2)
plot(time,Psol,'color',[0 0.5 0],'linewidth',3)
hold on
plot(time,smooth(Pmot,0.1),'r','linewidth',3)
set(gca,'fontsize',14)
grid on
legend('Available Solar Power','Required Motor Power')
ylabel('Power (W)')
xlabel('Time (hours)')
ylim([0,max(Psol.*1.5)])
xlim([0,max(time)])
saveas(gcf,[resDir '\' 'Sim - Energy and Power'],'fig')
saveas(gcf,[resDir '\' 'Sim - Energy and Power.png'],'png')
%print([resDir '\' 'Sim - Energy and Power'],'-dpng','-r300')

% Trajectory
figure
plot3(xto,yto,zto,'r','linewidth',2)
hold on
plot3(xch,ych,zch,'g','linewidth',2)
plot3(xcl,ycl,zcl,'k','linewidth',2)
plot3(xcr,ycr,zcr,'linewidth',2)
% plot3(xp,yp,zp,'xk','linewidth',2)
grid on
axis equal
zlim([0,max(zp).*1.2])
set(gca,'fontsize',14)
lgd = legend('Takeoff','Charging Cruise','Neutral Climb','Max Charge Cruise',...
    'Location','SouthOutside','Orientation','horizontal');
set(lgd,'FontSize',8)
xlabel('Downrange (ft.)')
ylabel('Crossrange (ft.)')
zlabel('Altitude (ft.)')
title(sprintf('SUAV Flight Trajectory\nLaunch: %s | %d Hour Flight | [%s,%s]',...
    Mission.dateStart,Mission.flightHours,...
    num2str(Mission.location(1)),num2str(Mission.location(2))))
saveas(gcf,[resDir '\' 'Sim - Flight Trajectory'],'fig')
saveas(gcf,[resDir '\' 'Sim - Flight Trajectory.png'],'png')
%print([resDir '\' 'Sim - Flight Trajectory'],'-dpng','-r300')

%% Contour Mapping Plots
% Excess Time
figure
hold on
contour(spanVec,mbatVec,tmet,[0,1],'b','linewidth',2)
plot(spanVec(spanLoc),mbatVec(batLoc),'go','linewidth',2)

pcolor(spanVec,mbatVec,texc_plot);
shading interp;
colormap hot;
colorbar
caxis([0,max(max(texc_plot))])
colormap(flipud(colormap));
set(gca,'FontSize',14)
xlabel('Wingspan (m)')
ylabel('Battery Mass (kg)')
title(sprintf('Excess Time (hours)\nLaunch: %s | %d Hour Flight | [%s,%s]',...
    Mission.dateStart,Mission.flightHours,...
    num2str(Mission.location(1)),num2str(Mission.location(2))))

hold on
contour(spanVec,mbatVec,tmet,[0,1],'b','linewidth',2)
plot(spanVec(spanLoc),mbatVec(batLoc),'go','linewidth',2)
legend('Feasible Space','Maximum Score Solution','Location','SouthEast')

saveas(gcf,[resDir '\' 'Sizing - Excess Time'],'fig')
print([resDir '\' 'Sizing - Excess Time'],'-dpng','-r300')

% Charge Margin
figure
hold on
contour(spanVec,mbatVec,tmet,[0,1],'b','linewidth',2)
plot(spanVec(spanLoc),mbatVec(batLoc),'go','linewidth',2)

pcolor(spanVec,mbatVec,tcm_plot);
colormap hot;
shading interp;
colorbar
caxis([0,max(max(tcm_plot))])
colormap(flipud(colormap));
set(gca,'FontSize',14)
xlabel('Wingspan (m)')
ylabel('Battery Mass (kg)')
title(sprintf('Charge Margin (hours)\nLaunch: %s | %d Hour Flight | [%s,%s]',...
    Mission.dateStart,Mission.flightHours,...
    num2str(Mission.location(1)),num2str(Mission.location(2))))

hold on
contour(spanVec,mbatVec,tmet,[0,1],'b','linewidth',2)
plot(spanVec(spanLoc),mbatVec(batLoc),'go','linewidth',2)
legend('Feasible Space','Maximum Score Solution','Location','Northwest')

saveas(gcf,[resDir '\' 'Sizing - Charge Margin'],'fig')
print([resDir '\' 'Sizing - Charge Margin'],'-dpng','-r300')

% Time Metric
figure
hold on
plot(spanVec(spanLoc),mbatVec(batLoc),'go','linewidth',2)

pcolor(spanVec,mbatVec,tmet)
shading interp
colormap hot;
colorbar
colormap(flipud(colormap));
set(gca,'FontSize',14)
xlabel('Wingspan (m)')
ylabel('Battery Mass (kg)')
title(sprintf('Energy Robustness Score (Time Metric)\nLaunch: %s | %d Hour Flight | [%s,%s]',...
    Mission.dateStart,Mission.flightHours,...
    num2str(Mission.location(1)),num2str(Mission.location(2))))

hold on
plot(spanVec(spanLoc),mbatVec(batLoc),'go','linewidth',2)
legend('Maximum Score Solution','Location','Northwest')

saveas(gcf,[resDir '\' 'Sizing - Time Metric'],'fig')
print([resDir '\' 'Sizing - Time Metric'],'-dpng','-r300')

%% Save Max Score Aircraft
save([resDir '\maxScoreAircraft.mat'],'maxScoreAircraft');
save([resDir '\mission.mat'],'Mission');