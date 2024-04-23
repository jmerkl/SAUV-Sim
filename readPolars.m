function [alfas,Cl,Cd,Cdp,Cm,Xcp] = readPolars(airfoil)
selectNew = false;
localDesig = 'Airfoils\';

matfile = dir([localDesig airfoil '.mat']);
if isempty(matfile) || selectNew
    polfiles = dir([localDesig '*.txt']);
    for i = 1:length(polfiles)
       filename = polfiles(i).name;
       if ~isempty(strfind(filename, airfoil))
          filenames{i} =  filename;
       end
    end
    fileLoc = 1;
    if length(filenames) > 1
        fileLoc = menu('Select Desired Polar File', filenames);
    end
    filename = filenames{fileLoc};
    
    fh = fopen([localDesig filename],'r');
    line = fgetl(fh);
    polFound = false;
    j = 1;
    while ischar(line) %end of the text file
        if polFound
            lineCopy = line;
            i = 1;
            while ~isempty(lineCopy)
               [val, lineCopy] =  strtok(lineCopy);
               polVec(i) = str2num(val);
               i = i + 1;
            end
            polMat(j,:) = polVec;
            j = j + 1;
        end
        if ~isempty(strfind(line,'------- -------- ---------'))
           polFound = true; 
        end
        line = fgetl(fh);
    end
    alfas = deg2rad(polMat(:,1));
    Cl = polMat(:,2);
    Cd = polMat(:,3);
    Cdp = polMat(:,4);
    Cm = polMat(:,5);
    Xcp = polMat(:,10);
    save([localDesig airfoil '.mat'],'alfas','Cl','Cd','Cdp','Cm','Xcp');
end
load([localDesig airfoil '.mat']);