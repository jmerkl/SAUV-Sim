%clear
clc
close all

%% Load Battery File
batteryData;
% Or xlsread & format the values from an Excel sheet

%% Individual SoC Current Plotting
mkdir(batteryName); clc;
colors = {'b','r','g','k','y','m','c'};
[r,c] = size(Vmat);
Itest = Imat;
for i = 1:r
   Imat(i,1:c) = Imat(i,1) .* ones(1,c);
   Qpmat(i,1:c) = (Qmat(i,1:c) ./ Qmat(i,c)) .* 100;
   
   Qvec = spline(1:length(Qpmat(i,:)),Qpmat(i,:),linspace(1,length(Qpmat),1e3));
   Vvec = spline(Qpmat(i,:),Vmat(i,:),Qvec);
   Vmati(i,:) = Vvec;
   Qpmati(i,:) = Qvec;
   Qmati(i,:) = Qvec .* Qmat(i,end);
   
   hold on
   plot((Qpmati(i,:)./100).*Qmati(i,end)./360000, Vmati(i,:),colors{i},'LineWidth',2)
   legendCell{i} = sprintf('%sA',num2str(Imat(i)));
end
grid on
set(gca,'FontSize',14)
legend(legendCell)
xlabel('Battery Capacity in A-hr')
ylabel('Cell Voltage')
title(sprintf('%s: State of Charge',batteryName))
print(sprintf('%s/%s SoC.png',batteryName,batteryName),'-dpng','-r300')
saveas(gca,sprintf('%s/%s SoC.fig',batteryName,batteryName))

%% Curve-to-Curve Interpolation across Current
Idis = linspace(min(min(Itest')),max(max((Itest'))),1e3);
Qmax = spline(Itest',Qmati(:,end),Idis);
[r,c] = size(Vmati);
for i = 1:c
    Vcont(i,:) = interp1(Qmati(:,end),Vmati(:,i),Qmax);
    Qcont(i,:) = interp1(Qmati(:,end),Qmati(:,i),Qmax);
end
Vcont = fliplr(Vcont');
figure
Qperc = linspace(0,100,c);
pcolor(Qperc, Idis, Vcont)
shading interp
colorbar
set(gca,'FontSize',14)
title('Cell SoC Map (Voltage)')
ylabel('Discharge Current in Amps')
xlabel('Battery Charge (in %)')
print(sprintf('%s/%s SoC Map.png',batteryName,batteryName),'-dpng','-r300')
saveas(gca,sprintf('%s/%s SoC Map.fig',batteryName,batteryName))

varBattName = batteryName;
varBattName(strfind(varBattName,' ')) = '';
save([batteryName '/' varBattName '.mat'],'Vcont','Idis','Qmax','Qcont');