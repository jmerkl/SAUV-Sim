function [Aircraft] = massEnergySizing(Aircraft, Mission)
%% Unpack Variables
% Battery Parameters
Vcell = Aircraft.Batteries.Vcell;
Ccell = Aircraft.Batteries.Ccell;
Cmin = Aircraft.Batteries.Cmin;
mcell = Aircraft.Batteries.BatteryCell.mcell;

% Solar Parameters
FFsol = Aircraft.Solar.FF;
Acell = Aircraft.Solar.SolarCell.Acell;
insolAve = (2/pi) .* max(Mission.Solar.Insol);
ksm = Aircraft.Solar.SolarCell.ksm;
kmppt = Aircraft.Solar.kmppt;
n_mod = Aircraft.Solar.n_mod;
n_mppt = Aircraft.Solar.n_mppt;

% Wing Parameters
b = Aircraft.Aero.Wing.span;
AR = Aircraft.Aero.Wing.AR;

% Mission Parameters
Vto = Mission.Vto;
Kto = Mission.Kto;
tn = Mission.tnight;
td = Mission.tday;
r = Mission.Dloiter ./ 2;

% Propulsion Parameters
Vmot = Aircraft.Propulsion.Vmot;
n_prop = Aircraft.Propulsion.n_prop;

% Mass Values
mbat = Aircraft.Mass.batteries;
maf = Aircraft.Mass.airframe;
mpl = Aircraft.Mass.payload;
mav = Aircraft.Mass.avionics;
msm = Aircraft.Mass.solar;
m_mppt = Aircraft.Mass.mppt;
mprop = Aircraft.Mass.propulsion;

%% Configure Battery Pack
if Mission.discrete
    numSer = ceil(Vmot ./ Vcell);
    numPar = floor((mbat ./ mcell) ./ numSer);
else
    numSer = (Vmot ./ Vcell);
    numPar = ((mbat ./ mcell) ./ numSer);
end

Vpack = numSer .* Vcell;
Cpack = numPar .* Ccell;
Epack = Vpack .* Cpack;

%% Initialize Total Mass
mo = mbat + maf + mpl + mprop + mav + msm + m_mppt;

%% Payload/Avionics Power
Paux = Aircraft.Electronics.Pav + Aircraft.Electronics.Ppl;

%% Generate Takeoff Standard Atmosphere
[Atmosphere] = calcAtmosphere(Mission.hto);
rhoTo = Atmosphere.rho;
go = Atmosphere.g;

%% Propulsion Mass Estimation (Takeoff Condition)
CLmax = Aircraft.Aero.Polars.CLmax;
CDmax = Aircraft.Aero.Polars.CDmax;
theta_to = Aircraft.Aero.Polars.alfaMax;
qto = 0.5 .* rhoTo .* (Vto ./ Kto).^2;
Sto = Aircraft.Aero.Wing.Swet;
Lto = qto .* Sto .* CLmax;
Dto = qto .* Sto .* CDmax;
Tto_x = (Dto.*cos(theta_to)) + (Lto.*sin(theta_to)) + (mo .* go .* sin(theta_to));
Tto_y = (mo .* go) - Lto.*cos(theta_to) + Dto.*sin(theta_to);
Tto = norm([Tto_x Tto_y]);
Pto = Tto .* Vto; %Takeoff power required - should be the maximum motor demand
Pmot_max = Pto ./ n_prop;

mprop = 0.0011 .* Pmot_max;

%% Generate Gliding Standard Atmosphere
[Atmosphere] = calcAtmosphere(Mission.hmax);
rhoGl = Atmosphere.rho;
go = Atmosphere.g;

%% Solar Climb/Glide Sizing
CDcr = Aircraft.Aero.Polars.CDcr;
CLcr = Aircraft.Aero.Polars.CLcr;
dEpack = Epack .* (1 - (Cmin./1e2));

Psol = @(S) insolAve .* S .* (1 ./ FFsol) .* n_mod .* n_mppt;
Pcr = @(S) (0.5 .* rhoGl .* CDcr .* S .* ((2 .* mo .* go) ./ (S .* rhoGl .* CLcr)).^(3/2))./n_prop;
fA = @(S) (tn ./ (td - (dEpack ./ (Psol(S) - Pcr(S)))));
fB = @(S) (1 ./ (mo .* go)) .* ((Psol(S) .* n_prop) - Pcr(S));

Svec = linspace(0, 2 .* Sto);
for i = 1:length(Svec)
    fAv(i) = fA(Svec(i));
    fBv(i) = fB(Svec(i));
end
dhdt = fBv ./ fAv;
dh = dhdt .* tn;
[~,loc] = min(abs(dh - Mission.hmax));
SsolBound = Svec(loc);

dfAdS = gradient(fAv,Svec);
dfBdS = gradient(fBv,Svec);
d2hdt2 = (fAv .* dfBdS) - (fBv .* dfAdS);
[~,loc] = min(abs(d2hdt2));
SsolUnbound = Svec(loc);

%% Choose Largest Required Area
%SwetVec = max([SsolBound SsolUnbound Sto]);
Swet = Aircraft.Aero.Wing.Swet;

%% Solar Mass Estimation
chord = (Swet ./ b) .* sqrt(FFsol);
if Mission.discrete
    numVer = floor(chord ./ sqrt(Acell));
    numHor = floor((b./2) ./ sqrt(Acell)) .* 2;
else
    numVer = (chord ./ sqrt(Acell));
    numHor = ((b./2) ./ sqrt(Acell)) .* 2;
end
numTot = numVer .* numHor;
Asol = Acell .* numTot;
Psm = Psol(Swet);
msm = ksm .* Asol;

%% MPPT Mass Estimation
m_mppt = kmppt .* Psm;

%% Airframe Mass Estimation
maf = (0.44./go) .* b.^3.1 .* AR.^-0.25;

%% Structure Mass Estimation
nmax = Aircraft.Struct.Wing.nmax;
mstr = (0.44 .* nmax./go) .* (Swet./maf).^3.1 .* (AR.^-0.25);

%% Total Mass Build-up
mtot = mbat + maf + mstr + mpl + mprop + mav + msm + m_mppt;

%% Set Up Minimum Power Cruise Equations
[~,rho] = calcAtmosphere((Mission.hmax + Mission.hcruise)./2);

alfas = Aircraft.Aero.Polars.alfas;
CD = Aircraft.Aero.Polars.CD_total;
CL = Aircraft.Aero.Polars.CL;

Vlev = @(gamma,CLlev) sqrt(2.*go.*mtot.*r) ./ ((CLlev.*rho.*r.*cos(gamma).*Swet).^2 - 4.*mtot.^2).^0.25;
phi = @(V) atan(V.^2 ./ (go.*r));
q = @(V) 0.5 .* rho .* V.^2;
Pmot = @(V,dhdt,CD) (q(V).*Swet.*CD.*V + mtot.*go.*dhdt) .* (1./n_prop) + Paux;
CLreq = @(V,gamma) (mtot .*go .* (1 + sin(gamma)) .* (1./cos(phi(V))))./(q(V).*Swet);

%% Calculate Excess Time
Pbat = dEpack./tn;
dhdt_gl = -(Mission.hmax - Mission.hcruise)./tn; %negative to represent descent

Vs = Vlev(0,CL);
Pnet = Pmot(Vs,dhdt_gl,CD)  - Pbat;
dPdV = gradient(Pnet,Vs);
[~,locV] = min(abs(dPdV));
Vcr = Vs(locV);
CLcr = CL(locV);
CDcr = CD(locV);
alfaCr = alfas(locV);

gamma_gl = asin(dhdt_gl ./ Vcr);
Vreq = Vlev(gamma_gl,CLcr);
Vgl = min([Vreq Vcr]);
gamma_gl = asin(dhdt_gl ./ Vgl);

CLgl = CLreq(Vgl,gamma_gl);
[~,glLoc] = min(abs(CL - CLgl));
CLgl = CL(glLoc);
CDgl = CD(glLoc);
alfaGl = alfas(glLoc);

dEdt = Pmot(Vgl,dhdt_gl,CDgl);
Enight = dEdt .* tn;
Epf = Epack - Enight;
Epf(Epf < 0) = 0;

texc = Epf ./ dEdt;

%% Calculate Charge Margin Time
if Epf < (Epack - dEpack)
    Eo = Epf; %Initial pack energy as the amount left at the end of the night
else
    Eo = Epack - dEpack; %Initial pack energy as the starting amount
end

Pncr = Psm - Pmot(Vcr,0,CDcr);
tch = (Epack - Eo) ./ Pncr;
tch(tch < 0) = 0;

dhdt_cl_min = (Mission.hmax - Mission.hcruise) ./ (td - tch);
dhdt_cl_max = Vcr .* sin(Aircraft.Aero.Polars.alfaMax);
dhdts = linspace(dhdt_cl_min,dhdt_cl_max);
gammas = asin(dhdts./Vcr);

CLcls = CLreq(Vcr,gammas);
for i = 1:length(CLcls)
    [~,clLoc] = min(abs(CL - CLcls(i)));
    CDcls(i) = CD(clLoc);
end
Vcls = Vlev(gammas,CLcls);
Pncl = Psm - Pmot(Vcls,dhdts,CDcls);
[~,minLoc] = min(abs(Pncl));

Vcl = Vcls(minLoc);
gamma_cl = gammas(minLoc);

CLcl = CLreq(Vcl,gamma_cl);
[~,clLoc] = min(abs(CL - CLcl));
CDcl = CD(clLoc);
alfaCl = alfas(clLoc);

dhdt_cl = Vcl .* sin(gamma_cl);
tcl = (Mission.hmax - Mission.hcruise) ./ dhdt_cl;

tcm = td - tch - tcl;
tcm(tcm < 0) = 0;

%% Repack Variables
% Aircraft Performance
Aircraft.Performance.chargeMargin = tcm;
Aircraft.Performance.excessTime = texc;
Aircraft.Performance.cruiseSpeed = Vcr;
Aircraft.Performance.climbSpeed = Vcl;
Aircraft.Performance.glideSpeed = Vgl;
Aircraft.Performance.cruiseAoA = rad2deg(alfaCr);
Aircraft.Performance.climbAoA = rad2deg(alfaCl);
Aircraft.Performance.glideAoA = rad2deg(alfaGl);
Aircraft.Performance.cruiseFPA = rad2deg(0);
Aircraft.Performance.climbFPA = rad2deg(gamma_cl);
Aircraft.Performance.glideFPA = rad2deg(gamma_gl);
Aircraft.Performance.cruisePitch = rad2deg(alfaCr);
Aircraft.Performance.climbPitch = rad2deg(alfaCl) + rad2deg(gamma_cl);
Aircraft.Performance.glidePitch = rad2deg(alfaGl) + rad2deg(gamma_gl);

% Battery Pack Parameters
Aircraft.Batteries.Epack = Epack;
Aircraft.Batteries.Vpack = Vpack;
Aircraft.Batteries.Cpack = Cpack;
Aircraft.Batteries.series = numSer;
Aircraft.Batteries.parallel = numPar;

% Solar Array Parameters
Aircraft.Solar.area = Asol;
Aircraft.Solar.vertical = numVer;
Aircraft.Solar.horizontal = numHor;

% Takeoff Wing Area
Aircraft.Aero.Wing.Swet = Swet;
Aircraft.Aero.Wing.chord = chord;
Aircraft.Aero.Wing.AR = AR;

% Mass Parameters
Aircraft.Mass.airframe = maf;
Aircraft.Mass.structure = mstr;
Aircraft.Mass.solar = msm;
Aircraft.Mass.mppt = m_mppt;
Aircraft.Mass.propulsion = mprop;
Aircraft.Mass.total = mtot;
end