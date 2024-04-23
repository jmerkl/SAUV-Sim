function [Atmosphere, rho] = calcAtmosphere(alt)
%% Sea Level Constants
Po = 101325;
To = 288.15;
Htro = 10972.8; %average height of the troposphere

% Air Specific Constants
R = 287.058; %J/(kg-K)
C1 = 1.458e-6; %Sutherland coefficient
S = 110.4; %Sutherland Temperature
gamma = 1.4; %ratio of specific heats for air

% Non-flat earth constants
ME = 5.972e24; %kg
RE = 6.371e6; %m
G = 6.674e-11; %Gravitational Constant

%% Round-Earth Correction
g = G .* (ME ./ (RE + alt).^2);

%% 1976 Atmosphere Model Calculation
% Invalid altitude measurements
if alt > 100000
    error('Units not in meters and/or altitude is in space.  No space.')
end
% Temperature
T = To - (6.5 .* alt./1000);
T(T < 216.65) = 216.65;
% Pressure
P11 = Po .* (1 - (0.0065 .* Htro ./ To)).^5.2561;
T11 = To - (6.5 .* Htro./1000);

P = P11 .* exp((-g./(R .* T11)) .* (alt - Htro));
P(alt < Htro) = Po .* (1 - (0.0065 .* alt ./ To)).^5.2561;

% Density
rho = P ./ (R .* T);
% Kinematic Viscosity
u = (C1 .* T.^(3/2)) ./ (T + S); %Sutherland's Law
kvs = u ./ rho; %kinematic viscosity
% Speed of Sound
a = sqrt(R .* gamma .* T);

%% Package Outputs
Atmosphere.temp = T;
Atmosphere.press = P;
Atmosphere.rho = rho;
Atmosphere.kvs = kvs;
Atmosphere.a = a;
Atmosphere.g = g;
end