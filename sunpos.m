function PSAplus = sunpos(TS, Lat, Long)
%
% Inputs: 
%       - TS                : Time stamp (Nx6) vector in the order Year,
%                           Month, Day, hour, minute second, as one would
%                           obtaine with datevec command
%       - Lat, Long         : Latitude and Longitude, in degrees
%
%
% Outputs: 
%       - PSA   : Structure with calculated variables,
%               + ELong:            Ecliptic Longitude, in radians
%               + ELong_deg:        Ecliptic Longitude, in degrees
%               + EObl:             Ecliptic Obliquity, in radians
%               + EObl_deg:         Ecliptic Obliquity, in degrees
%               + RightAscension:   Right Ascension, in radians
%               + Declination:      Declination, in radians
%               + HourAngle:        Hour angle, in radians
%               + Azimuth:          Solar Azimuth angle, from North, in degrees
%               + Azimuth_rad:      Solar Azimuth angle, from North, in radians
%               + Zenith:           Solar Zenith angle, in degrees
%               + Zenith_rad:       Solar Zenith angle, in radians
%               + SunVec:           [Nx3] sun vector in {east, north, zenith} 
%                                   coordinate system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate Julian Date and Elapsed Julian Days
[EJD, DecimalHours] = ElapsedJulianDays(TS);

% PSA Algorithm

%% Ecliptic Coordinates
dOmega =  2.267127827e+00 - 9.300339267e-04*EJD;
dMeanLongitude = 4.895036035e+00 + 1.720279602e-02*EJD; 
dMeanAnomaly = 6.239468336e+00 + 1.720200135e-02*EJD;

% Ecliptic longitude
PSAplus.ELong = dMeanLongitude + 3.338320972e-02*sin( dMeanAnomaly ) ...
    + 3.497596876e-04 * sin( 2*dMeanAnomaly ) - 1.544353226e-04 ...
    - 8.689729360e-06*sin( dOmega );
PSAplus.ELong_deg = PSAplus.ELong * 180/pi;

% Ecliptic Obliquity
PSAplus.EObl = 4.090904909e-01 - 6.213605399e-09*EJD + ...
    4.418094944e-05*cos(dOmega);
PSAplus.EObl_deg = PSAplus.EObl * 180/pi;


%% Celestial coordinates
% Right ascension and declination
dY1 = cos( PSAplus.EObl ) .* sin( PSAplus.ELong );
dX1 = cos( PSAplus.ELong );
PSAplus.RightAscension = atan2( dY1, dX1);
ind = find(PSAplus.RightAscension < 0);
PSAplus.RightAscension(ind) = PSAplus.RightAscension(ind) + 2*pi;
PSAplus.Declination = asin( sin( PSAplus.EObl ) .* sin( PSAplus.ELong ) );


%% Topocentric coordinates
% Greenwich and Local sidereal time, Hour angle
dGreenwichMeanSiderealTime = 6.697096103e+00 + 6.570984737e-02*EJD + DecimalHours;
dLocalMeanSiderealTime = ( dGreenwichMeanSiderealTime*15+Long )*pi/180;
PSAplus.HourAngle = dLocalMeanSiderealTime - PSAplus.RightAscension;

% Local coordinates
PSAplus.Zenith_rad = (acos( cos(Lat*pi/180)*cos( PSAplus.HourAngle ).*cos( PSAplus.Declination ) + ...
    sin( PSAplus.Declination )*sin(Lat*pi/180)));
dY2 = -sin(PSAplus.HourAngle);
dX2 = tan( PSAplus.Declination) * cos( Lat*pi/180 ) - sin(Lat*pi/180)*cos(PSAplus.HourAngle);
PSAplus.Azimuth_rad = atan2(dY2, dX2);
ind = find(PSAplus.Azimuth_rad < 0);
PSAplus.Azimuth_rad(ind) = PSAplus.Azimuth_rad(ind) + 2*pi;

% Parallax correction
PSAplus.Zenith_rad = PSAplus.Zenith_rad + 6371.01/149597870.7 * sin(PSAplus.Zenith_rad);


% convert to degrees
PSAplus.Zenith = rad2deg(PSAplus.Zenith_rad);
PSAplus.Azimuth = rad2deg(PSAplus.Azimuth_rad);


PSAplus.SunVec = [sin(PSAplus.Azimuth_rad).*sin(PSAplus.Zenith_rad), ...
    cos(PSAplus.Azimuth_rad).*sin(PSAplus.Zenith_rad), cos(PSAplus.Zenith_rad)];

end


function [EJD, DecimalHours] = ElapsedJulianDays(TS)
    Year = TS(:,1); Month = TS(:,2); Day = TS(:,3);
    Hour = TS(:,4); Minute = TS(:,5); Second = TS(:,6);

    DecimalHours = Hour + Minute/60 + Second/3600;

    JulianDate = fix( (1461*( Year +4800+ fix(( Month-14)/12) ))/4 ) + ...
        fix( (367 * (Month-2-12*fix((Month-14)/12)))/12 ) - ...
        fix( 3*( fix((Year+4900 + fix((Month-14)/12))/100))/4 ) + ...
        Day - 32075.5 + DecimalHours/24;

    % Elapsed Julian Days, since 2000 January 1 at 1200h
    EJD = JulianDate - 2451545.0;
end