clear
clc
close all

%% Maximum Power Point Tracker
MPPT.MPPTname = 'Common';
MPPT.mpptEff = @(boost) -1e-4 .* (boost).^2 - 0.004.*boost + 0.9941;

save(['Model Params\MPPTs\' MPPT.MPPTname '.mat']);
fprintf('%s Parameters Saved...\n',MPPT.MPPTname)
pause(0.1)
clear

%% Solar Cells
SolarCell.cellName = 'Sunpower C60';
SolarCell.n_sol = 0.226; %cell efficiency
SolarCell.Acell = 0.0155; %square meters
SolarCell.ksm = 0.59; %solar cell mass/area in kg/m^2 (includes EVA encapsulation)
SolarCell.Isc = 6.32; %short circuit current for a single cell
SolarCell.Io = 1e-10; %saturation current as determined by recombination
SolarCell.Voc = 0.68; %Cell open-circuit voltage
SolarCell.kTq = 0.02586; %Thermal voltage @ 300K
SolarCell.dVdT = -0.3055; % %/K
SolarCell.dIdT = 0.0455; %  %/K
SolarCell.dPdT = -0.3910; % %/K
SolarCell.Ttest = 25; %Tested temperature in Celsius
SolarCell.Tsat = 80;
SolarCell.n = 0.7; %Ideality factor

save(['Model Params\Solar Cells\' SolarCell.cellName '.mat']);
fprintf('%s Parameters Saved...\n',SolarCell.cellName)
pause(0.1)
clear

%% Battery Cells
% Battery cell information
BatteryCell.batteryName = 'Panasonic 18650b';
BatteryCell.socProfile = 'Model Params\Battery Cell Models\Panasonic 18650b\Panasonic18650b.mat'; %Location of the SoC profile
BatteryCell.mcell = 0.0475;
BatteryCell.Rcell = 1.43e-4; %Ohms/cell w/ tabbed connections

% Battery operation limits
BatteryCell.Vmax = 4.15;
BatteryCell.Vmin = 2.55;
BatteryCell.Imax = 20;
BatteryCell.Imin = -4;
BatteryCell.Tmax = 50;

% Load SoC Parameters
load(BatteryCell.socProfile);
BatteryCell.SoC.Idis = Idis;
BatteryCell.SoC.Qmax = Qmax;
BatteryCell.SoC.Qcont = Qcont;
BatteryCell.SoC.Vcont = Vcont;

save(['Model Params\Batteries\' BatteryCell.batteryName '.mat']);
fprintf('%s Parameters Saved...\n', BatteryCell.batteryName)
pause(0.1)
clear

%% Update Model List
paramsFolder = 'Model Params';
folders = dir(paramsFolder);
dirFlags = [folders.isdir];
subFolders = folders(dirFlags);
for i = 3:length(subFolders)
   folderNames{i-2} = subFolders(i).name;
end

fh = fopen([paramsFolder '\' 'Model_List.txt'],'wt');
fprintf(fh,'################# LIST OF ALL STORED MODEL PARAMETERS #################\n\n');
for i = 1:length(folderNames)
    headString = ['--------------' folderNames{i} '--------------\n\n'];
    fprintf(fh,headString);
    matFiles = dir([paramsFolder '\' folderNames{i} '\' '*.mat']);
    for j = 1:length(matFiles)
        matFile = matFiles(j).name;
        matFileDir = [paramsFolder '\' folderNames{i} '\' matFile];
        fprintf(fh,[matFile ':\n']);
        varNames = who('-file',matFileDir);
        vars = load(matFileDir,varNames{1});
        varFields = fieldnames(vars.(varNames{1}));
        vars = vars.(varNames{1});
        for k = 1:length(varFields)
           var = varFields{k};
           val = vars.(varFields{k});
           if strcmpi('double',class(val))
               val = num2str(val);
           elseif iscell(val)
               tempVal = [];
               for ct = 1:length(val)
                   if ct < length(val)
                       tempVal = [tempVal val{ct} ';'];
                   else
                       tempVal = [tempVal val{ct}];
                   end
               end
               val = tempVal;
           end
           if strcmpi(class(val),'function_handle'); val = func2str(val); end
           if ~isstruct(val)
                fprintf(fh,['     ' var ': ' val '\n']);
           end
        end
        fprintf(fh,'\n\n');
    end
end
fclose(fh);
fprintf('\nModel List Updated!\n\n')

clear