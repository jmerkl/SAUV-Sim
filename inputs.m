% Aircraft
Aircraft.Performance.staticMargin = -0.25;

Aircraft.Aero.Wing.e = 0.92; %Oswald efficiency factor
Aircraft.Struct.Fuselage.df = 0.5 * (12 / 39.6); %Required Fuselage diameter (ft. -> m)
Aircraft.Struct.Fuselage.lambdap = 2; %length of payload section as a ratio of fuselage diameter
Aircraft.Struct.Fuselage.lambdaf = 14; %fuselage fineness ratio
Aircraft.Aero.tail2wing = 0.3; %initial area ratio of the tail to the wing
Aircraft.Aero.Wing.airfoil = 'SF-15'; %Folder designation of airfoil polars
Aircraft.Aero.Tail.airfoil = 'NACA 0008';

% Structure
Aircraft.Struct.Wing.nmax = 7; %maximum normal g-loading

% Solar
Aircraft.Solar.solarCell = 'Sunpower C60';
Aircraft.Solar.n_mod = 0.23;
Aircraft.Solar.n_mppt = 0.95;
Aircraft.Solar.FF = 1.2; %Solar Fill Factor

% Batteries
Aircraft.Batteries.batteryCell = 'Panasonic 18650b';
Aircraft.Batteries.Vcell = 3.7; %Average voltage
Aircraft.Batteries.Ccell = 3.4 .* 3600; %Average capacity (A-s)
Aircraft.Batteries.Cmin = 10; %Minimum Charge Percentage (as %)
Aircraft.Batteries.Ccl = 95; %Charge to start climbing (as %)

% Electronics
Aircraft.Propulsion.n_prop = 0.6; %propulsion system efficiency
Aircraft.Solar.mppt = 'Common';
Aircraft.Solar.kmppt = 0.422e-3; %mass constant kg/W
Aircraft.Propulsion.Vmot = 22.2; %Motor Voltage
Aircraft.Electronics.Pav = 3; %W
Aircraft.Electronics.Ppl = 0.5; %W

% On-board masses
Aircraft.Mass.avionics = 0.5; %mass of avionics in kg
Aircraft.Mass.payload = 1;
Aircraft.Mass.airframe = 0;
Aircraft.Mass.solar = 0;
Aircraft.Mass.mppt = 0;
Aircraft.Mass.propulsion = 0;