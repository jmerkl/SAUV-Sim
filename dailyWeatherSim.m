function [InsolShade,T,C,ncc,t] = dailyWeatherSim(Insol,CC,CCtime,Tamb,Tsat,time)
%% Cloud Cover Calculation
t = time .* 60; %convert to minutes
% Interpolate Cloud Cover Points
%CCint = spline(CCtime,CC,linspace(0,length(Insol)));
%CCs = CCint((start*60):60:(finish*60)) ./ 100;
CCs = 0 .* ones(1,length(Insol));

% Calculate Shading Efficiency
ncc = 1 - 0.75.*CCs.^3.4; %scienceDirect model
InsolShade = Insol .* ncc;

%% Heating Loss Calculation
To = Tamb;
Tf = Tsat;
C = -log(Tf - To);
KT = (log(Tf - 1.2*To) - C) ./ time(end); %KT = 0.0001214;
T = Tf - exp(-((KT .* time) + C));

%% Factor in Cloud Cover to Heating
%ncShade = sqrt(ncc);
%T = T .* ncShade;
end