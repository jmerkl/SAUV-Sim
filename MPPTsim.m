function [Pmpp, Vmpp, Impp] = MPPTsim(solarArray, Irr, Isc, Voc, eff, Acell, series, parallel)
load(solarArray);
Area = series .* parallel .* Acell;
Prr = Irr .* Area .* eff;
IL = Prr ./ (series .* Voc);
if IL > Isc
    IL = Isc; %light generated current cannot be greater than the short circuit current
end

M = parallel; %# of cells in parallel
Vt = series * [0:(Voc ./ 100):Voc]; %total voltage from the module

i = 1;
cond = true;
while cond && i <= length(Vt)
    It(i) = (M .* IL) - (M .* Io .* (exp((Vt(i)) ./ (series .* kTq .* n)) - 1));
    P(i) = It(i) .* Vt(i);
    if P(i) < 0
       cond = false;
    end
    i = i + 1;
end

Impp = It(P == max(P));
Vmpp = Vt(P == max(P));
Pmpp = Impp .* Vmpp;

%% Data Export and Graphing
% figure
% subplot(2, 2, 1)
% set(gca,'FontSize',14)
% plot(Vt(1:i-2), P(1:i-2),'LineWidth',2)
% xlabel('Module Voltage in Volts')
% ylabel('Module Power in Watts')
% title('Module P-V Curve')
% grid on
% 
% subplot(2, 2, 3)
% set(gca,'FontSize',14)
% plot(Vt(1:i-2), It(1:i-2),'LineWidth',2)
% xlabel('Module Voltage in Volts')
% ylabel('Module Current in Amps')
% title('Module I-V Curve')
% axis([min(Vt), max(Vt), 0, max(It) * 1.1])
% grid on
% 
% subplot(2, 2, [2,4])
% [ax, p1, p2] = plotyy(Vt(1:i-2), P(1:i-2), Vt(1:i-2), It(1:i-2), 'plot');
% set(p1, 'LineWidth', 2)
% set(p2, 'LineWidth', 2)
% grid on
% title('Module MPP Curve')
% ylabel(ax(1),'Module Power in Watts')
% ylabel(ax(2),'Module Current in Amps')
% xlabel('Module Voltage in Volts')
% end