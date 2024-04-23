function [Vcur,Qcur,Qspent,effPack,lowestV,Vpsoc,Qpsoc] = batteryModel(packName, Pm, Psol, t, Qprev)
%% Load Variables
packFile = ['Model Params\Packs\' packName '.mat'];
load(packFile); %loads in vehicle config selected pack
load(['Battery Cell Models\' SoCname '\' SoCname '.mat']);

%% Establish Current SoC Line
Qfull = nPar .* Qcell;
Vnom = nSer .* Vcell;
if isnan(Qprev)
    Qprev = Qfull; %Full pack
end

Pnet = Pm - Psol; %Positive power leaves the pack, negative charges it
Inet = Pnet ./ (Vnom - Rpack.*(Pnet./Vnom));
Icell = mean(Inet(~isnan(Inet)) ./ nPar);
ind = find(Idis > Icell); ind = ind(1);
Cmax = Qmax(ind)./360000; Vsoc = Vcont(ind,:);
if Cmax > Qcell
    Cmax = Qcell;
end
Csoc = (linspace(0,1,length(Vsoc))).*Cmax;
Vpsoc = Vsoc .* nSer;
Qpsoc = Csoc .* nPar;

%% Calculate Dissipative Power Loss
Ploss = Inet.^2 .* Rpack;
effPack = 1 - (Ploss ./ Pnet);

%% Coulomb Count Current to Determine new SoC Position
dt = median(gradient(t));
Qspent = cumsum(Inet ./ effPack).*dt ./ 3600;
Qnet = Qspent(end); %Euler-style summation integration, A-hr
Qcur = Qprev - Qnet;
[~,loc] = min(abs(Qpsoc - Qcur));
Vcur = Vpsoc(loc);
Qperc = 100 .* (Qcur ./ Qfull); %Percent of Charge Left

%% Maximum Imbalance Calc
cellVs = (Vcur ./ nSer) .* ones(1,nSer);
Rmod = Rpack ./ nSer;
for i = 1:length(cellVs)
   celldIs = (Inet.^2 .* i.*Rmod)./cellVs(i);
   celldQs(i) = sum(celldIs(i).*dt); 
end
minQ = (Qprev./nPar) - min(celldQs)./3600;
maxQ = (Qprev./nPar) - max(celldQs)./3600;
[~,minLoc] = min(abs((Qpsoc./nPar)-minQ));
[~,maxLoc] = min(abs((Qpsoc./nPar)-maxQ));
minV = (Vpsoc(minLoc)./nSer);
maxV = (Vpsoc(maxLoc)./nSer);
dVcell = maxV - minV;
lowestV = maxV;
if lowestV < Vmin
    %error('BMS TRIP - UNDERVOLTAGE FROM IMBALANCE')
end
end