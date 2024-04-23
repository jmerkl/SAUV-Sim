function [Aircraft] = aeroPolarCalc(Aircraft, Atmosphere)
%% Initialize Sizing of Non-Wing Surfaces
[Aircraft] = initialSizing(Aircraft);

%% Unpack parameters (for easier readability)
Vinf = Aircraft.Performance.cruiseSpeed;

AR_w = Aircraft.Aero.Wing.AR;
AR_th = Aircraft.Aero.Tail.Horizontal.AR;
AR_tv = Aircraft.Aero.Tail.Vertical.AR;

b_w = Aircraft.Aero.Wing.span;
b_th = Aircraft.Aero.Tail.Vertical.span;
b_tv = Aircraft.Aero.Tail.Horizontal.span;

e_w = Aircraft.Aero.Wing.e;
e_th = Aircraft.Aero.Wing.e; %Assuming same oswald for tail as wing
e_tv = Aircraft.Aero.Wing.e;

Swf = Aircraft.Aero.Fuselage.Swf;
L_f = Aircraft.Struct.Fuselage.Lf;

t2w = Aircraft.Aero.tail2wing;

%% Wing 2-D to 3-D Calculation
[Aircraft.Aero.Wing.Polars] = liftingLineCalc(Aircraft.Aero.Wing.airfoil, Vinf, b_w, AR_w, e_w);
CDw = Aircraft.Aero.Wing.Polars.CD;
CLw = Aircraft.Aero.Wing.Polars.CL;

%% Horizontal Stabilizer Drag Build-up
[Aircraft.Aero.Tail.Horizontal.Polars] = liftingLineCalc(Aircraft.Aero.Tail.airfoil, Vinf, b_th, AR_th, e_th);

% Scale polars to wing reference area
Aircraft.Aero.Tail.Horizontal.Polars.CDo = Aircraft.Aero.Tail.Horizontal.Polars.CDo .* t2w;
Aircraft.Aero.Tail.Horizontal.Polars.CD = Aircraft.Aero.Tail.Horizontal.Polars.CDo...
    + (Aircraft.Aero.Tail.Horizontal.Polars.k .* (Aircraft.Aero.Tail.Horizontal.Polars.CL.* t2w).^2);

if length(Aircraft.Aero.Wing.Polars.alfas) < length(Aircraft.Aero.Tail.Horizontal.Polars.alfas)
    for i = 1:length(Aircraft.Aero.Wing.Polars.alfas)
       [~,minLoc] = min(abs(Aircraft.Aero.Wing.Polars.alfas(i) - Aircraft.Aero.Tail.Horizontal.Polars.alfas));
       CDth(i,1) = Aircraft.Aero.Tail.Horizontal.Polars.CD(i);
       CLth(i,1) = Aircraft.Aero.Tail.Horizontal.Polars.CL(i);
    end
    alfas = Aircraft.Aero.Wing.Polars.alfas;
    Aircraft.Aero.Tail.Horizontal.Polars.alfas = alfas;
else
    CDw = [];
    CLw = [];
    for i = 1:length(Aircraft.Aero.Tail.Horizontal.Polars.alfas)
       [~,minLoc] = min(abs(Aircraft.Aero.Wing.Polars.alfas(i) - Aircraft.Aero.Tail.Horizontal.Polars.alfas));
       CDw(i) = Aircraft.Aero.Wing.Polars.CD(i);
       CLw(i) = Aircraft.Aero.Wing.Polars.CL(i);
       CDth(i) = Aircraft.Aero.Tail.Horizontal.Polars.CD(i);
       CLth(i) = Aircraft.Aero.Tail.Horizontal.Polars.CL(i);
    end
    alfas = Aircraft.Aero.Tail.Horizontal.Polars.alfas;
    Aircraft.Aero.Wing.Polars.alfas = alfas;
end

Aircraft.Aero.Wing.Polars.CL = CLw';
Aircraft.Aero.Wing.Polars.CD = CDw';
Aircraft.Aero.Tail.Horizontal.Polars.CLw = CLth';
Aircraft.Aero.Tail.Horizontal.Polars.CDw = CDth';

%% Vertical Stabilizer Drag Build-up
[Aircraft.Aero.Tail.Vertical.Polars] = liftingLineCalc(Aircraft.Aero.Tail.airfoil, Vinf, b_tv, AR_tv, e_tv);
Aircraft.Aero.Tail.Vertical.Polars.CDo = Aircraft.Aero.Tail.Vertical.Polars.CDo .* (Aircraft.Aero.Tail.Vertical.Swet./Aircraft.Aero.Wing.Swet);
CDtv = Aircraft.Aero.Tail.Vertical.Polars.CDo;

%% Fuselage Drag Build-up
Re = Vinf .* L_f ./ Atmosphere.kvs;
Cdf = (0.074.*Re.^(-0.2));
CDf = Cdf .* (Swf ./ Aircraft.Aero.Wing.Swet); %Scale to wing reference area for build-up

Aircraft.Aero.Fuselage.CDf = CDf;

%% Total Drag Polar Estimation
CD = CDw + CDf + CDth + CDtv;
Aircraft.Aero.Polars.CD_total = CD';
Aircraft.Aero.Polars.CL = CLw';

%% Cruise and Max Lift Condition
DLe = CD ./ (CLw.^(3./2));

[~,loc] = min(DLe);   [~,maxLoc] = max(CLw);
alfaCr = alfas(loc);  alfaMax = alfas(maxLoc);
CLcr = CLw(loc);       CLmax = CLw(maxLoc);
CDcr = CD(loc);       CDmax = CD(maxLoc);

%% Parasitic Drag
[~,parLoc] = min(abs(alfas));
CDo = CD(parLoc);

%% Package Aero Results
Aircraft.Aero.Polars.alfas = alfas;

Aircraft.Aero.Polars.alfaCr = alfaCr;
Aircraft.Aero.Polars.CLcr = CLcr;
Aircraft.Aero.Polars.CDcr = CDcr;
Aircraft.Aero.Wing.Polars.Cmcr = Aircraft.Aero.Wing.Polars.Cm(loc);

Aircraft.Aero.Polars.alfaMax = alfaMax;
Aircraft.Aero.Polars.CLmax = CLmax;
Aircraft.Aero.Polars.CDmax = CDmax;
Aircraft.Aero.Polars.CDo = CDo;

end

function [Polars] = liftingLineCalc(airfoil, Vinf, b, AR, e)
%% Read in 2-D Airfoil Polars
[alfas,Cl,Cd,Cdp,Cm,Xcp] = readPolars(airfoil);

%% Coefficient Derivatives
[~,endLoc] = max(Cl.^(3/2) ./ Cd);
alfaEnd = alfas(endLoc);
Cla = mode(gradient(Cl(1:endLoc),alfas(1:endLoc)));

%% Calculate Induced AoA for each Wing AoA (for downwash)
LLon = false;

if LLon
    n = 20;
    s = b./2;
    theta0 = linspace(0,pi,n);
    y0 = -s.*cos(theta0);
    y = linspace(-s,s,length(theta0));
    c=interp1([-s,0,s],[1,1,1],y0); %chord distribution
    alpha0 = zeros(n,1);

    for k = 1:length(alfas)
        alpha = alfas(k).*ones(n,1);
        Ano = alpha(2:n-1)-alpha0(2:n-1);

        % Create Fourier Constants
        for i = 2:n-1
            for j = 2:n-1
                A(i-1,j-1)=4*s*sin((j-1)*theta0(i))/Cla/c(i)+(j-1)*sin((j-1)*theta0(i))/sin(theta0(i));
            end
        end
        An = A\Ano;

        % Calculate Vortex Strengths
        for i = 2:n-1
           Gamma(i)=2*b*Vinf*sum(An(:).*sin((1:n-2)*theta0(i))');
        end
        Gamma(n) = 0;
        Gamma(1) = 0;

        % Calculate Downwash Angles
        for i=2:n-1,
            alphai(i)=sum((1:n-2)'.*An.*sin((1:n-2)'*theta0(i))./sin(theta0(i)));
        end
        alphai(1)=sum((1:n-2)'.*An.*(1:n-2)');
        alphai(n)=sum((1:n-2)'.*An.*(1:n-2)');

        Polars.alfai(k,:) = alphai;
    end
    Polars.yspan = y;
else
   Polars.alfai = ((2.*Cla)./(pi.*e.*AR)).*alfas; 
end

%% 2-D -> 3-D PLL Calculation of Lift Curve
dCldA = gradient(Cl,alfas);
dA = gradient(alfas);
dCLdA = dCldA * (AR ./ (AR + 2));
CL = cumtrapz(dCLdA .* dA) + Cl(1);

%% Correct Polars with CL
k = 1 ./ (pi .* e .* AR);
CDo = Cd(alfas == 0);
CDi = k .* CL.^2;
CD = Cd + CDi;

%% Package Output
Polars.alfas = alfas;
Polars.k = k;
Polars.Cl = Cl;
Polars.Cd = Cd;
Polars.Cdp = Cdp;
Polars.Cm = Cm;
Polars.Xcp = Xcp;
Polars.CL = CL;
Polars.CDi = CDi;
Polars.CDo = CDo;
Polars.CD = CD;
end

function [Aircraft] = initialSizing(Aircraft)
%% First Loop Initialize Tail Sizing
if ~isfield(Aircraft.Aero.Tail, 'Horizontal')
    ARw = Aircraft.Aero.Wing.AR;
    ARt = 2 ./ ((2).^(1 - ARw./(ARw+2)) .* ((ARw + 2)./ARw) - 1);
    Aircraft.Aero.Tail.Horizontal.AR = ARt;
    Aircraft.Aero.Tail.Horizontal.Swet = Aircraft.Aero.Wing.Swet .* Aircraft.Aero.tail2wing;
    Aircraft.Aero.Tail.Horizontal.span = sqrt(ARt.*Aircraft.Aero.Tail.Horizontal.Swet);
    Aircraft.Aero.Tail.Horizontal.chord = Aircraft.Aero.Tail.Horizontal.Swet ./ Aircraft.Aero.Tail.Horizontal.span;

    Aircraft.Aero.Tail.Vertical.chord = Aircraft.Aero.Tail.Horizontal.chord;
    Aircraft.Aero.Tail.Vertical.span = Aircraft.Aero.Tail.Horizontal.span ./ 2;
    Aircraft.Aero.Tail.Vertical.Swet = Aircraft.Aero.Tail.Vertical.chord .* Aircraft.Aero.Tail.Vertical.span;
    Aircraft.Aero.Tail.Vertical.AR = Aircraft.Aero.Tail.Vertical.span.^2 ./ Aircraft.Aero.Tail.Vertical.Swet;
end

%% Initial Fuselage & Payload Section Sizing
d_f = Aircraft.Struct.Fuselage.df;
lambda_f = Aircraft.Struct.Fuselage.lambdaf;
lambda_p = Aircraft.Struct.Fuselage.lambdap;

L_f = lambda_f .* d_f;
L_p = lambda_p .* d_f;
Swf = pi .* d_f .* L_f .* ((1 - (2./lambda_f)).^(2/3)) .* (1 + (1 ./ lambda_f.^2));

Aircraft.Aero.Fuselage.Swf = Swf;
Aircraft.Struct.Fuselage.Lf = L_f;
Aircraft.Struct.Fuselage.Lp = L_p;

end