function [Aircraft] = calcSuavInertiaCG(Aircraft)
%% Initialize Variables
Lf = Aircraft.Struct.Fuselage.Lf;
Lp = Aircraft.Struct.Fuselage.Lp;
Aircraft.Aero.Wing.thickness = Aircraft.Aero.Wing.chord .* 0.15;
Aircraft.Aero.Tail.Horizontal.thickness = Aircraft.Aero.Tail.Horizontal.chord .* 0.08;
Aircraft.Aero.Tail.Vertical.thickness = Aircraft.Aero.Tail.Vertical.chord .* 0.08;

%% Calculate Mass Locations
% X-Position
Aircraft.Struct.Xpos.Tail.Horizontal = Lf;
Aircraft.Struct.Xpos.Tail.Vertical = Lf;
Aircraft.Struct.Xpos.Fuselage = Lf ./ 2;
Aircraft.Struct.Xpos.Batteries = Lp;
Aircraft.Struct.Xpos.Wing = Lp + 0.5.*Aircraft.Aero.Wing.chord;
Aircraft.Struct.Xpos.Solar = Aircraft.Struct.Xpos.Wing;
Aircraft.Struct.Xpos.Payload = Lp ./ 2;
Aircraft.Struct.Xpos.Propulsion = 0;

% Y-Position (Assuming a left-right symmetrical aircraft)
Aircraft.Struct.Ypos.Tail.Horizontal = 0;
Aircraft.Struct.Ypos.Tail.Vertical = 0;
Aircraft.Struct.Ypos.Fuselage = 0;
Aircraft.Struct.Ypos.Batteries = 0;
Aircraft.Struct.Ypos.Wing = 0;
Aircraft.Struct.Ypos.Solar = 0;
Aircraft.Struct.Ypos.Payload = 0;
Aircraft.Struct.Ypos.Propulsion = 0;

% Z-Position (vert./horiz. stabilizers)
ztv = (Aircraft.Struct.Fuselage.df ./ 2) + (Aircraft.Aero.Tail.Vertical.span ./ 2);
zth = (Aircraft.Struct.Fuselage.df ./ 2) + (Aircraft.Aero.Tail.Vertical.span ./ 2) + (Aircraft.Aero.Tail.Horizontal.thickness./2);

Aircraft.Struct.Zpos.Tail.Horizontal = zth;
Aircraft.Struct.Zpos.Tail.Vertical = ztv;
Aircraft.Struct.Zpos.Fuselage = 0;
Aircraft.Struct.Zpos.Batteries = 0;
Aircraft.Struct.Zpos.Wing = 0;
Aircraft.Struct.Zpos.Solar = 0;
Aircraft.Struct.Zpos.Payload = 0;
Aircraft.Struct.Zpos.Propulsion = 0;

%% Calculate CG Location
maf = Aircraft.Mass.airframe;
mbat = Aircraft.Mass.batteries;
mpl = Aircraft.Mass.payload;
mav = Aircraft.Mass.avionics;
msm = Aircraft.Mass.solar;
m_mppt = Aircraft.Mass.mppt;
mprop = Aircraft.Mass.propulsion;

Xpos = [Aircraft.Struct.Xpos.Tail.Horizontal, Aircraft.Struct.Xpos.Tail.Vertical,...
        Aircraft.Struct.Xpos.Fuselage, Aircraft.Struct.Xpos.Batteries,...
        Aircraft.Struct.Xpos.Wing, Aircraft.Struct.Xpos.Solar,...
        Aircraft.Struct.Xpos.Payload, Aircraft.Struct.Xpos.Propulsion];
    
Ypos = [Aircraft.Struct.Ypos.Tail.Horizontal, Aircraft.Struct.Ypos.Tail.Vertical,...
        Aircraft.Struct.Ypos.Fuselage, Aircraft.Struct.Ypos.Batteries,...
        Aircraft.Struct.Ypos.Wing, Aircraft.Struct.Ypos.Solar,...
        Aircraft.Struct.Ypos.Payload, Aircraft.Struct.Ypos.Propulsion];

Zpos = [Aircraft.Struct.Zpos.Tail.Horizontal, Aircraft.Struct.Zpos.Tail.Vertical,...
        Aircraft.Struct.Zpos.Fuselage, Aircraft.Struct.Zpos.Batteries,...
        Aircraft.Struct.Zpos.Wing, Aircraft.Struct.Zpos.Solar,...
        Aircraft.Struct.Zpos.Payload, Aircraft.Struct.Zpos.Propulsion];

Stot = Aircraft.Aero.Wing.Swet + Aircraft.Aero.Fuselage.Swf + Aircraft.Aero.Tail.Horizontal.Swet + Aircraft.Aero.Tail.Vertical.Swet;
mw = maf .* (Aircraft.Aero.Wing.Swet ./ Stot);
mth = maf .* (Aircraft.Aero.Tail.Horizontal.Swet ./ Stot);
mtv = maf .* (Aircraft.Aero.Tail.Vertical.Swet ./ Stot);
mf = maf .* (Aircraft.Aero.Fuselage.Swf ./ Stot);

Mvec = [mth, mtv, mf, mbat, mw, msm, (mpl + mav + m_mppt), mprop];

Xcg = sum(Xpos .* Mvec) ./ sum(Mvec); %X-axis Center of gravity from the nose of the aircraft
Xg = Xcg - (Aircraft.Struct.Xpos.Wing - Aircraft.Aero.Wing.chord./4); %Wing AC to CG distance
Xt = Aircraft.Struct.Xpos.Tail.Horizontal - Xcg; %Tail to CG distance

Ycg = sum(Ypos .* Mvec) ./ sum(Mvec); %Y-Axis Center of gravity from the fuselage centerline

Zcg = sum(Zpos .* Mvec) ./ sum(Mvec); %Z-axis Center of gravity from the fuselage centerline

%% Set Up Inertia Tensor Equations

% Approximate Wing Surface Inertia Equations (Rectangular Prism)
Ixxs = @(m,t,b) (1./12).*m.*(t.^2 + b.^2);
Iyys = @(m,t,c) (1./12).*m.*(t.^2 + c.^2);
Izzs = @(m,c,b) (1./12).*m.*(c.^2 + b.^2);

% Approximate Fuselage Inertia Equations Equations (Solid Rod)
Ixxr = @(m,r) 0.5.*m.*r.^2;
Iyyr = @(m,r,L) (1./12).*m.*(3.*r.^2 + L.^2);
Izzr = Iyyr;

%% Calculate Inertia Tensor
% Calculate Wing Inertia Values
mwi = mw + msm; %wing inertia mass includes wing airframe and solar module weight
Ixxw = Ixxs(mwi,Aircraft.Aero.Wing.thickness,Aircraft.Aero.Wing.span);
Iyyw = Iyys(mwi,Aircraft.Aero.Wing.thickness,Aircraft.Aero.Wing.chord);
Izzw = Izzs(mwi,Aircraft.Aero.Wing.chord,Aircraft.Aero.Wing.span);

% Calculate Horizontal Stabilizer Inertia Values
Ixxh = Ixxs(mth,Aircraft.Aero.Tail.Horizontal.thickness,Aircraft.Aero.Tail.Horizontal.span);
Iyyh = Iyys(mth,Aircraft.Aero.Tail.Horizontal.thickness,Aircraft.Aero.Tail.Horizontal.chord);
Izzh = Izzs(mth,Aircraft.Aero.Tail.Horizontal.chord,Aircraft.Aero.Tail.Horizontal.span);

% Calculate Vertical Stabilizer Inertia Values
Ixxv = Ixxs(mtv,Aircraft.Aero.Tail.Vertical.thickness,Aircraft.Aero.Tail.Vertical.span);
Iyyv = Izzs(mtv,Aircraft.Aero.Tail.Vertical.thickness,Aircraft.Aero.Tail.Vertical.chord); %Iyy and Izz flip due to vert. vs horiz.
Izzv = Iyys(mtv,Aircraft.Aero.Tail.Vertical.chord,Aircraft.Aero.Tail.Vertical.span);

% Calculate Fuselage Inertia Values
mfi = [mf  mprop  mbat  (mpl + mav + m_mppt)];
Ixxf = Ixxr(sum(mfi),Aircraft.Struct.Fuselage.df./2);
Iyyf = Iyyr(sum(mfi),Aircraft.Struct.Fuselage.df./2,Aircraft.Struct.Fuselage.Lf);
Izzf = Izzr(sum(mfi),Aircraft.Struct.Fuselage.df./2,Aircraft.Struct.Fuselage.Lf);

% Sum Polar Inertia Values
miVec = [sum(mfi) mw mth mtv];

IxxVec = [Ixxf Ixxw Ixxh Ixxv];
XfVec = [Aircraft.Struct.Xpos.Fuselage, Aircraft.Struct.Xpos.Propulsion,...
         Aircraft.Struct.Xpos.Batteries, Aircraft.Struct.Xpos.Payload];
Xf = sum(XfVec .* mfi) ./ sum(mfi);
XxVec = [Xf, Aircraft.Struct.Xpos.Wing, Aircraft.Struct.Xpos.Tail.Horizontal, Aircraft.Struct.Xpos.Tail.Vertical];
Ixx = sum(IxxVec +  miVec.*abs(XxVec - Xcg).^2); %parallel axis theorem

IyyVec = [Iyyf Iyyw Iyyh Iyyv];
YyVec = [Aircraft.Struct.Ypos.Fuselage, Aircraft.Struct.Ypos.Wing,...
         Aircraft.Struct.Ypos.Tail.Horizontal, Aircraft.Struct.Ypos.Tail.Vertical];
Iyy = sum(IyyVec + miVec.*abs(YyVec - Ycg));

IzzVec = [Izzf Izzw Izzh Izzv];
ZzVec = [Aircraft.Struct.Zpos.Fuselage, Aircraft.Struct.Zpos.Wing,...
         Aircraft.Struct.Zpos.Tail.Horizontal, Aircraft.Struct.Zpos.Tail.Vertical];
ZzVec = sum(IzzVec + miVec.*abs(ZzVec - Zcg));
Izz = sum(IzzVec + miVec.*abs(ZzVec - Zcg));

% Calculate products of inertia (Ixz, Izx, etc.)
Ixy = sum(miVec .* abs(XxVec - Xcg) .* abs(YyVec - Ycg));
Ixz = sum(miVec .* abs(XxVec - Xcg) .* abs(ZzVec - Zcg));
Iyz = sum(miVec .* abs(YyVec - Ycg) .* abs(ZzVec - Zcg));
Izx = Ixz;
Iyx = Ixy;
Izy = Iyz;

% Assemble Inertia Tensor
I = [Ixx -Ixy -Ixz
    -Iyx  Iyy -Iyz
    -Izx -Izy  Izz];

%% Package Output Variables
Aircraft.Mass.CG.Xcg = Xcg;
Aircraft.Mass.CG.Ycg = Ycg;
Aircraft.Mass.CG.Zcg = Zcg;
Aircraft.Mass.InertiaTensor = I;
Aircraft.Stability.Xg = Xg;
Aircraft.Stability.Xt = Xt;
end