clear
clc
close all

% Simulation
dt = 15;

% Solar Information
nSol = 0.23;
aSol = 1.24;

% Location
Mission.dateStart = '2023/6/15 0:00:00';
Mission.location = [30.22, -92.03]; %lat/long position
Mission.altitude = 0 .* (39.6 / 12);
Mission.flightHours = 72; %number of consecutive hours of flight
Mission.nightFactor = 1.4; %increase the duration of "night" where Ps < Pmot

time = 0:dt:(Mission.flightHours * 3600);
[Insol, Az, El] = solarCurveSim(Mission.location(1),Mission.location(2),Mission.altitude,Mission.dateStart,time);
Psol = nSol .* aSol .* Insol .* sin(deg2rad(El));


subplot(3,1,1)
plot(time./3600,El,'linewidth',2)
grid on
set(gca,'fontsize',14)
ylabel('Elevation Angle (degrees)')
subplot(3,1,2)
plot(time./3600,Insol,'linewidth',2)
grid on
set(gca,'fontsize',14)
ylabel('Local Insolation (W/m^2)')
xlabel('Time (hours)')
subplot(3,1,3)
plot(time./3600,Psol,'linewidth',2)
grid on
set(gca,'fontsize',14)
ylabel('Solar Power (W)')
xlabel('Time (hours)')

