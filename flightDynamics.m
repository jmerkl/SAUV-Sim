% clear
% clc
% close all

%% Simulate Full Flight
% Unpack Variables
Aircraft = maxScoreAircraft;
Vpack = Aircraft.Batteries.Vpack;
Cpack = Aircraft.Batteries.Epack ./ Vpack;
C = Mission.chargeStart .* Cpack;
r = Mission.Dloiter;

% Pre Calculate Raw Solar
Insol = Insol - min(Insol);
Psol = Insol .* (Aircraft.Solar.n_mod .* Aircraft.Solar.n_mppt) .* Aircraft.Solar.area;

% Initialize Metrics
dEdt = [0 0];
depletedBattery = false;
TOAchieved = false;

% Initialize Position
z = Mission.hto;
y = 0;
x = 0;

for k = 1:length(t)
    Cp = (C(k) ./ Cpack) .* 1e2;
    Ps = Psol(k);

    % Modes of Operation
    dayTime = Ps > 0;
    nightTime = Ps <= 0;
    charging = dEdt(k) > 0;

    if Cp <= 0
        depletedBattery = true;
    end

    if ~depletedBattery
        takeoffMode(k) = (z(k) < Mission.hcruise) && ~TOAchieved;
        chargingCruiseMode(k) = ((Cp <= Aircraft.Batteries.Ccl) || (z(k) >= Mission.hmax)) && charging && dayTime;
        neutralClimbMode(k) = (Cp > Aircraft.Batteries.Ccl) && (z(k) < Mission.hmax) && dayTime;
        nightGlideMode(k) = (~charging) && (z(k) > Mission.hcruise) || nightTime;
        
        if ~takeoffMode(k) && ~neutralClimbMode(k) && ~nightGlideMode(k)
           chargingCruiseMode(k) = true; 
        end
    end

    % Generate Atmosphere
    [Atmosphere] = calcAtmosphere(z(k));
    rho = Atmosphere.rho;
    Ta = Atmosphere.temp;
    g = Atmosphere.g;

    % Extract Values
    mtot = Aircraft.Mass.total;
    Swet = Aircraft.Aero.Wing.Swet;
    n_prop = Aircraft.Propulsion.n_prop;

    % Simulate Conditions
    if depletedBattery
        CL = Aircraft.Aero.Polars.CLmax;
        CD = Aircraft.Aero.Polars.CDmax;
        alfa(k) = Aircraft.Aero.Polars.alfaMax;
        gamma(k) = -atan(CD./CL);

        theta(k) = alfa(k) + gamma(k);
        phi(k) = 0;
        psi(k) = 0;

        Vac(k) = 0;
        q(k) = 0.5 .* rho .* Vac(k).^2;
        Dac(k) = q(k) .* Swet .* CD;
        Lac(k) = q(k) .* Swet .* CL;

        Pmot(k) = 0;

    elseif takeoffMode(k) %Takeoff Condition
        Vac(k) = Aircraft.Performance.climbSpeed;     
        q(k) = 0.5 .* rho .* Vac(k).^2;
        CLreq = (mtot.*g) ./ (q(k) .* Swet);
        [~,toLocCL] = min(abs(Aircraft.Aero.Polars.CL - CLreq));
        
        alfa(k) = Aircraft.Aero.Polars.alfas(toLocCL);
        CL = Aircraft.Aero.Polars.CL(toLocCL);
        CD = Aircraft.Aero.Polars.CD_total(toLocCL);
        
        gamma(k) = deg2rad(Aircraft.Performance.climbFPA);
        
        Lac(k) = q(k) .* Swet .* CL;
        Dac(k) = q(k) .* Swet .* CD;
        Tto_x = (Dac(k).*cos(gamma(k))) + (Lac(k).*sin(gamma(k))) + (mtot.*g.*sin(gamma(k)));
        Tto_y = (mtot.*g) - Lac(k).*cos(gamma(k)) + Dac(k).*sin(gamma(k));
        T(k) = norm([Tto_x Tto_y]);

        Pmot(k) = T(k) .* Vac(k) ./ n_prop;

        % Calculate Body Angles
        theta(k) = alfa(k) + gamma(k);
        psi(k) = 0;
        phi(k) = 0;
        
    elseif chargingCruiseMode(k) % Charging Cruise
        CL = Aircraft.Aero.Polars.CLcr;
        CD = Aircraft.Aero.Polars.CDcr;
        alfa(k) = Aircraft.Aero.Polars.alfaCr;
        gamma(k) = 0;

        Vac(k) = sqrt(2 .* g .* mtot .* r) ./ ((CL.^2 .* rho.^2 .* r.^2 .* Swet.^2) - 4.*mtot.^2).^0.25;
        phi(k) = atan(Vac(k).^2 ./ (g .* r));                

        q(k) = 0.5 .* rho .* Vac(k).^2;
        Lac(k) = q(k) .* Swet .* CL;
        Dac(k) = q(k) .* Swet .* CD;
        T(k) = Dac(k) + mtot.*g.*sin(gamma(k));
        Pmot(k) = T(k) .* Vac(k) ./ n_prop;

        psidot = Vac(k) ./ r;

        psi(k) = psidot.*dt + psi(k-1);
        theta(k) = alfa(k) + gamma(k);
        
        TOAchieved = true;

    elseif neutralClimbMode(k) %Neutral Climb
        
        CL = Aircraft.Aero.Polars.CLcr;
        CD = Aircraft.Aero.Polars.CDcr;
        alfa(k) = Aircraft.Aero.Polars.alfaCr;

        maxGamma = asin(T(1) ./ (mtot.*g));
        gammas = linspace(0,maxGamma);
        Vs = sqrt(2 .* g .* mtot .* r) ./ ((CL.^2 .* rho.^2 .* r.^2 .* cos(gammas).^2 .* Swet.^2) - 4.*mtot.^2).^0.25;
        qcl = 0.5 .* rho .* Vs.^2;
        Preq = (((qcl .* Swet .* CD) + mtot.*g.*sin(gammas)) .* Vs) ./ n_prop;
        dP = Ps - Preq;
        [~,loc] = min(abs(dP));

        Vac(k) = Vs(loc);
        Pmot(k) = Preq(loc);

        gamma(k) = gammas(loc);
        theta(k) = gamma(k) + alfa(k);
        psidot = Vac(k) ./ r;
        psi(k) = psidot.*dt + psi(k-1);
        phi(k) = atan(Vac(k).^2 ./ (g .* r));

        q(k) = 0.5 .* rho .* Vac(k).^2;
        Lac(k) = q(k) .* Swet .* CL;
        Dac(k) = q(k) .* Swet .* CD;
        T(k) = Dac(k) + + mtot.*g.*sin(gamma(k));
        
        TOAchieved = true;

     elseif nightGlideMode(k) %Night powered glide
         
        CL = Aircraft.Aero.Polars.CLcr;
        CD = Aircraft.Aero.Polars.CDcr;
        alfa(k) = Aircraft.Aero.Polars.alfaCr;
        gamma(k) = deg2rad(Aircraft.Performance.glideFPA);
        theta(k) = gamma(k) + alfa(k);

        Vac(k) = sqrt(2 .* g .* mtot .* r .* cos(gamma(k))) ./...
            ((CL.^2 .* rho.^2 .* r.^2 .* Swet.^2) - 4.*mtot.^2).^0.25;

        phi(k) = atan(Vac(k).^2) ./ (g .* r);
        psidot = Vac(k) ./ r;
        psi(k) = psidot.*dt + psi(k-1);

        q(k) = 0.5 .* rho .* Vac(k).^2;
        Lac(k) = q(k) .* Swet .* CL;
        Dac(k) = q(k) .* Swet .* CD;
        T(k) = Dac(k) + mtot.*g.*sin(gamma(k));
        if T(k) < 0; T(k) = 0; end
        Pmot(k) = T(k) .* Vac(k) ./ n_prop;

        TOAchieved = true;
    end

    Pac(k) = Pmot(k) + Aircraft.Electronics.Ppl + Aircraft.Electronics.Pav;
    Pnet(k) = Ps - Pac(k);
    dEdt(k+1) = Pnet(k);
    Imot = Pnet(k) ./ Vpack;
    dCmot = Imot .* dt;
    C(k+1) = C(k) + dCmot;

    C(C > Cpack) = Cpack;
    C(C < 0) = 0;

    u(k) = Vac(k) .* cos(gamma(k)) .* cos(psi(k));
    v(k) = Vac(k) .* sin(psi(k)) .* cos(phi(k));
    w(k) = Vac(k) .* sin(gamma(k));

    x(k+1) = u(k).*dt + x(k);
    y(k+1) = v(k).*dt + y(k);
    z(k+1) = w(k).*dt + z(k);

    z(z < 0) = 0;
    
    frac = k ./ length(t);
    waitbar(frac)
end
