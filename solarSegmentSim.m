function [Psol, Isol, Vsol, time] = solarSegmentSim(Aircraft,solarArray, mpptFile, Vpack, dateStart, flightHours, location, maxInsol, Ta)
% Solar Array & MPPT Specifications
load(solarArray);
%load(packFile);
load(mpptFile);

%% Load Information
% Position
lat = location(1);
lon = location(2);
altitude = location(3);

% Weather and Environments
Tamb = Ta - 273.15;
CCtime = 0:flightHours;
CC = 0 .* ones(1,length(CCtime));

%% Simulate Solar Curves
[InsolPure, time, ~, ~] = solarCurveSim(lat,lon,altitude,dateStart,flightHours,maxInsol);
[Insol,Tsol,~,~,t] = dailyWeatherSim(InsolPure,CC,CCtime,Tamb,Tsat,time);

%% Simulate Maximum Power Point Tracking
for i = 1:length(Insol)
    for j = 1:length(seriesVec)
        dT = abs(Tsol(i) - Ttest);
        Isc_T = Isc .* (1 - (dIdT .* dT ./ 100));
        Voc_T = Voc .* (1 - (dVdT .* dT ./ 100));
        [Pmpp(i,j), Vmpp(i,j), Impp(i,j)] = MPPTsim(solarArray, Insol(i), Isc_T, Voc_T, n_sol, Acell, seriesVec(j), 1);
    end
end
Vsol = sum(Vmpp');
n_mppt = mpptEff(Vpack ./ Vsol); %MPPT boost ratio efficiency
Psol = sum(Pmpp') .* n_mppt;
Isol = Psol ./ Vpack;
end