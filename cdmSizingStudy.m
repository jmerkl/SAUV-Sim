clear
clc
close all

%% INPUTS
% Dynamic Variables
spanVec = 2.5:0.1:7; %Wingspan ranges (min/max, m)
mbatVec = [0 6]; %Battery pack mass ranges (min/max, kg)
Aircraft.Aero.Wing.AR = 18;

% Simulation Variables
discrete = true; %keeps battery/solar cells as whole units if true
convErr = 0.01; %solution convergence error
dt = 1; %discrete time increment for flight simulation (s)
dtSol = 10; %discrete time increment for solar simulation (s)

% Mission Profile
Mission.dateStart = '2023-6-21 6:00:00';
Mission.location = [30.22, -92.03]; %lat/long position
Mission.flightHours = 72; %number of consecutive hours of flight
Mission.nightFactor = 1.1; %increase the duration of "night" where Ps < Pmot
Mission.chargeStart = 0.3; %initial starting battery charge (%, as decimal)

% Flight Profile
Mission.hto = 49 * (12 / 39.6); %altitude above sea level (ft. -> meters)
Mission.hcruise = 149 * (12 / 39.6); %cruise altitude (ft. -> meters)
Mission.hmax = 400 * (12 / 39.6); %maximum allowable altitude (ft. -> meters)
Mission.Dloiter = 1000 * (12 / 39.6); %loiter circle diameter (ft. -> meters)
Mission.Vto = 6; %hand-launched takeoff speed (m/s)
Mission.Kto = 1.2; %takeoff velocity factor (Vto = Kto * Vstall)

%% Direct Preset Load Files
% Load Inputs
inputs;
%motorFile = ['Model Params\Motors\' motor '.mat'];
packFile = ['Model Params\Batteries\' Aircraft.Batteries.batteryCell '.mat'];
solarFile = ['Model Params\Solar Cells\' Aircraft.Solar.solarCell '.mat'];
mpptFile = ['Model Params\MPPTs\' Aircraft.Solar.mppt '.mat'];
load(packFile); load(solarFile); load(mpptFile);

Aircraft.Batteries.Meta.packFile = packFile;
Aircraft.Solar.Meta.solarFile = solarFile;
Aircraft.Solar.Meta.mpptFile = mpptFile;

Aircraft.Solar.SolarCell = SolarCell;
Aircraft.Batteries.BatteryCell = BatteryCell;

%% Generate Atmosphere
[Atmosphere] = calcAtmosphere((Mission.hcruise + Mission.hmax)./2);
g = Atmosphere.g;

%% Generate Time Vectors
t = 0:dt:(Mission.flightHours * 3600);  %Flight time
tSol = 0:dtSol:t(end);                  %Solar calculation time

%% Generate Sunlight Profiles
[Insol, Az, El] = solarCurveSim(Mission.location(1),Mission.location(2),...
    Mission.hcruise,Mission.dateStart,tSol);
% Interpolate to flight time vector
Insol = interp1(tSol, Insol, t);
Az = interp1(tSol, Az, t);
El = interp1(tSol, El, t);

% Assign to Mission variable
Mission.Solar.Insol = Insol;
Mission.Solar.time = t;
Mission.Solar.Position.Az = Az;
Mission.Solar.Position.El = El;

% Define Average Day/Night Cycle
locDay = find(abs(gradient(Insol)) > 0);
locNight = find(abs(gradient(Insol)) <= 0);
nDays = Mission.flightHours ./ 24;
Mission.tday = (dt .* length(locDay)) ./ nDays;
Mission.tnight = ((dt .* length(locNight)) ./ nDays) .* Mission.nightFactor;

%% Calculate Battery Pack Module Increments
Mission.discrete = discrete;
mcellVec = mbatVec(1):Aircraft.Batteries.BatteryCell.mcell:mbatVec(end);
if Mission.discrete
    % Realistically increment by modules rather than individual cells
    numSer = ceil(Aircraft.Propulsion.Vmot ./ Aircraft.Batteries.Vcell);
    [numPar,batLoc] = unique(floor((mcellVec ./ Aircraft.Batteries.BatteryCell.mcell) ./ numSer));
    mbatVec = mcellVec(batLoc);
else
    mbatVec = mcellVec;
end

%% Iterate Sizing Loop
initAircraft = Aircraft; %Save initial aircraft parameters
ct = 0;
mtot = zeros(length(mbatVec), length(spanVec));
tcm = zeros(length(mbatVec), length(spanVec));
texc = zeros(length(mbatVec), length(spanVec));
for i = 1:length(mbatVec)
    for j = 1:length(spanVec)
        tStart = tic;
        Aircraft = initAircraft; %Reset to start aircraft parameters
        
        Aircraft.Mass.batteries = mbatVec(i);
        Aircraft.Aero.Wing.span = spanVec(j);
        
        %% Loop Energy-Based Mass Buildup
        convCond = false;
        mtotVec = [0 0];
        while ~convCond
            mtotVec(1) = mtotVec(2);
            Aircraft = sizeStaticStability(Aircraft, Mission);
            Aircraft = massEnergySizing(Aircraft, Mission);
            mtotVec(2) = Aircraft.Mass.total;
            err = abs((mtotVec(1) - mtotVec(2)) ./ mtotVec(2));
            
            if err < convErr
                convCond = true;
            end
        end
        mtot(i,j) = mtotVec(2);
        sizedAircraft(i,j) = Aircraft;

        %% Extract Time Metrics
        tcm(i,j) = sizedAircraft(i,j).Performance.chargeMargin;
        texc(i,j) = sizedAircraft(i,j).Performance.excessTime;
        
        %% Display Progress
        ct = ct + 1;
        frac = ct ./ (length(spanVec) .* length(mbatVec));
        waitbar(frac)
        
        timesPerLoop(ct) = toc(tStart);
        if mod(ct,5) == 0 || ct == 1
            timePerLoop = mean(timesPerLoop);
            clc
            fprintf('Estimated time remaining to complete: %s seconds\n', num2str(round(timePerLoop .* ((length(spanVec).*length(mbatVec)) - ct))))
        end
    end
end

%% Calculate Scoring Metric and Max-Score-Aircraft
texc_plot = texc ./ 3600;
tcm_plot = tcm ./ 3600;
mtot_n = min(min(mtot)) ./ mtot;
texc_n = texc ./ max(max(texc));
tcm_n = tcm ./ max(max(tcm));
ftmet = (max(max(tcm)) ./ max(max(texc)));
tmet_n = tcm_n .* texc_n.^(ftmet) .* mtot_n;
tmet = tmet_n ./ max(max(tmet_n));

tmetMax = max(max(tmet));
[batLoc,spanLoc] = find(tmet == tmetMax);

maxScoreAircraft = sizedAircraft(batLoc,spanLoc);
if isempty(batLoc) && isempty(spanLoc)
   maxScoreAircraft = sizedAircraft(end,end); 
end

%% Flight Simulation
flightDynamics;

%% Plot Results
plotPerformanceFigures;