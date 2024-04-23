function [ Aircraft ] = sizeStaticStability( Aircraft, Mission )
%% Generate Atmosphere
[Atmosphere] = calcAtmosphere((Mission.hmax + Mission.hcruise)./2);

%% Initialize Surface Geometry
% Wing Geometry
Aircraft.Aero.Wing.Swet = Aircraft.Aero.Wing.span.^2 ./ Aircraft.Aero.Wing.AR;
Aircraft.Aero.Wing.chord = Aircraft.Aero.Wing.span ./ Aircraft.Aero.Wing.AR;

% First-Loop Initialize Tail Surfaces Cruise-Speed
if ~isfield(Aircraft.Performance, 'cruiseSpeed')
    Aircraft.Performance.cruiseSpeed = Mission.Vto;
    Aircraft = initTailSurfaceSizing(Aircraft);
end

%% Set Up Sizing Loop
convCond = false;
convErr = 0.1;
Xcg = [0 0];

while ~convCond
    %% Call lifting-line aero polar calculations
    [Aircraft] = aeroPolarCalc(Aircraft, Atmosphere);

    %% Calculate Center of Gravity and Inertia Tensor
    [Aircraft] = calcSuavInertiaCG(Aircraft);

    %% Re-Size Tail Around Static Margin and Coefficient Targets
    [Aircraft] = calcLongStaticStability(Aircraft, Atmosphere);
    [Aircraft] = calcLatDirectStability(Aircraft, Mission, Atmosphere);

    %% Re-Calculate CG
    [Aircraft] = calcSuavInertiaCG(Aircraft);
    
    %% Convergence Criteria
    Xcg(1) = Xcg(2);
    Xcg(2) = Aircraft.Stability.Xcg;
    err = abs((Xcg(2) - Xcg(1)) ./ Xcg(2));
    
    if err < convErr
       convCond = true; 
    end
end
end

function [Aircraft] = calcLongStaticStability(Aircraft, Atmosphere)
%% Set Up Variables
rho = Atmosphere.rho;
Vac = Aircraft.Performance.cruiseSpeed;

Sw = Aircraft.Aero.Wing.Swet;
St = Aircraft.Aero.Tail.Horizontal.Swet;
ARw = Aircraft.Aero.Wing.AR;

Xcg = Aircraft.Mass.CG.Xcg;
Xg = Aircraft.Stability.Xg;
Xt = Aircraft.Stability.Xt;

c = Aircraft.Aero.Wing.chord;
SMd = Aircraft.Performance.staticMargin; %desired static margin

alfas = Aircraft.Aero.Wing.Polars.alfas;
CLw = Aircraft.Aero.Wing.Polars.CL;
CLt = Aircraft.Aero.Tail.Horizontal.Polars.CLw;

%% Moment Rate Balancing
Svec = linspace(0.1*St,Sw,1e3);
validLoc = find(alfas > 0 & alfas < max(alfas)/4);
dCLdAw = mean(gradient(CLw(validLoc),alfas(validLoc)));
dCLdAt = mean(gradient(CLt(validLoc),alfas(validLoc)));

q = 0.5.* rho .* Vac.^2;
dLwdA = q .* Sw .* dCLdAw;
dLtdA = q .* Svec .* dCLdAt;
SMv = ((dLwdA .* (Xg./c)) - (dLtdA .* (Xt./c))) ./ (dLwdA + dLtdA);

[~,locSM] = min(abs(SMv - SMd));
SM = SMv(locSM);
St = Svec(locSM);

%% Moment Coefficient Calculation
Xnp = Xcg - (SM .* c);
Vh = (St.*Xt) ./ (Sw.*c);
dEpsdA = (2.*dCLdAw)./(pi.*Aircraft.Aero.Wing.AR);
dMdA = (dLwdA.*Xg) - (dLtdA(locSM).*(1-dEpsdA).*Xt);
dCMdA = dMdA ./ (q.*Sw.*c);

%% Re-Calculate Tail Size
ARt = 1 ./ ((dCLdAw./dCLdAt).^(1 - ARw./(ARw+2)) .* ((ARw + 2)./ARw) - 1);
Aircraft.Aero.tail2wing = St ./ Sw;
Aircraft.Aero.Tail.Horizontal.Swet = St;
Aircraft.Aero.Tail.Horizontal.AR = ARt;
Aircraft.Aero.Tail.Horizontal.span = sqrt(Aircraft.Aero.Tail.Horizontal.AR .* St);
Aircraft.Aero.Tail.Horizontal.chord = St ./ Aircraft.Aero.Tail.Horizontal.span;

%% Assign Outputs
Aircraft.Aero.Tail.Horizontal.Vh = Vh;
Aircraft.Stability.Xnp = Xnp;
Aircraft.Stability.Xcg = Xcg;
Aircraft.Stability.SM = SM;
Aircraft.Stability.Cma = dCMdA;

end

function [Aircraft] = calcLatDirectStability(Aircraft, Mission, Atmosphere)
%% Set Up Variables
Sw = Aircraft.Aero.Wing.Swet;
b = Aircraft.Aero.Wing.span;
c = Aircraft.Aero.Wing.chord;
Lv = Aircraft.Stability.Xt;
Xg = Aircraft.Stability.Xg;
r = Mission.Dloiter ./ 2;
Vac = Aircraft.Performance.cruiseSpeed;

CLv = Aircraft.Aero.Tail.Vertical.Polars.CL;
CLw = Aircraft.Aero.Wing.Polars.CL;
alfas = Aircraft.Aero.Tail.Vertical.Polars.alfas;

%% Approximate Required Vertical Tail Volume
bv = sin(max(alfas)) .* Aircraft.Struct.Fuselage.Lf;
Cv = Aircraft.Aero.Tail.Horizontal.chord;
Sv = Cv .* bv;
ARv = bv.^2 ./ Sv;

Vv = (Sv.*Lv)./(Sw.*b);

%% Approximate Required Wing Dihedral
g = Atmosphere.g;
gamDih = 2.*atan(Vac.^2 ./ (r .* g));

%% Approximate Directional Coefficients
validLoc = find(alfas > 0 & alfas < max(alfas)/4);
% Cn_beta (sideslip directional stability)
kf = 1 - sqrt(1 ./ Aircraft.Struct.Fuselage.lambdaf); %approximation of fuselage instability contribution
dCLdAv = mode(gradient(CLv(validLoc),alfas(validLoc)));
dSigdBeta = 0; %assuming no significant sidewash contribution
Cnb = kf .* dCLdAv .* (1 - dSigdBeta) .* Vv;

% Cl_beta (sideslip roll stability)
dCLdAw = mode(gradient(CLw(validLoc),Aircraft.Aero.Wing.Polars.alfas(validLoc)));
Clbw = - sin(gamDih).*dCLdAw.*(2.*Xg./b);
zvt = bv ./ 2;
Clbv = - (zvt ./ Lv) .* Cnb;
Clb = Clbw + Clbv;

% Cl_phi (roll angle stability)
Clr = -(1/2) .* (b/c) .* dCLdAw .* (1-cos(gamDih));

%% Set Variables
Aircraft.Aero.Tail.Vertical.Vv = Vv;
Aircraft.Aero.Tail.Vertical.Swet = Sv;
Aircraft.Aero.Tail.Vertical.chord = Cv;
Aircraft.Aero.Tail.Vertical.span = bv;
Aircraft.Aero.Tail.Vertical.AR = ARv;

Aircraft.Aero.Wing.dihedral = gamDih;

Aircraft.Stability.Cnb = Cnb;
Aircraft.Stability.Clb = Clb;
Aircraft.Stability.Clr = Clr;
end

function [Aircraft] = initTailSurfaceSizing(Aircraft)
%% Horizontal Tail Volume Initial Sizing
Aircraft.Aero.Tail.Horizontal.Swet = Aircraft.Aero.Wing.Swet .* Aircraft.Aero.tail2wing;
Aircraft.Aero.Tail.Horizontal.chord = Aircraft.Aero.Wing.chord ./ 2;
Aircraft.Aero.Tail.Horizontal.span = Aircraft.Aero.Tail.Horizontal.Swet ./ Aircraft.Aero.Tail.Horizontal.chord;
Aircraft.Aero.Tail.Horizontal.AR = Aircraft.Aero.Tail.Horizontal.span.^2 ./ Aircraft.Aero.Tail.Horizontal.Swet;

%% Vertical Stabilizer Initial Sizing
Aircraft.Aero.Tail.Vertical.chord = Aircraft.Aero.Tail.Horizontal.chord;
Aircraft.Aero.Tail.Vertical.span = Aircraft.Aero.Tail.Horizontal.span ./ 2;
Aircraft.Aero.Tail.Vertical.Swet = Aircraft.Aero.Tail.Vertical.span .* Aircraft.Aero.Tail.Vertical.chord;
Aircraft.Aero.Tail.Vertical.AR = Aircraft.Aero.Tail.Vertical.span.^2 ./ Aircraft.Aero.Tail.Vertical.Swet;
end