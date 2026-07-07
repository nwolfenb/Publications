% Na_vs_f_v3.m
clear all; close all; clc

%% Path
addpath('..\..\IcyRF')

%% Defaults
colors2 = colororder;
ylim1 = [0 50];
ylim2 = [0 150];
fontsize = 8;

%% Ice Core
% Table 2, Grimm et al. (2015)
data = readcell('Cole-Cole.xlsx');
core_name = data(2:end,1);
depth = round(cell2mat(data(2:end,2)));

% Isolate Shared Data
ind = find(ismember(core_name, 'Vostok (Meteoric)') | ismember(core_name, 'Vostok (Accreted 1)'));

core_name = core_name(ind);
T0 = cell2mat(data(2:end,3));
T0 = T0(ind);
eps_inf = cell2mat(data(2:end,4));
eps_inf = eps_inf(ind);
delta_eps1 = cell2mat(data(2:end,5));
delta_eps1 = delta_eps1(ind);
fr1 = cell2mat(data(2:end,6));
fr1 = fr1(ind);
alpha1 = cell2mat(data(2:end,7));
alpha1 = alpha1(ind);
delta_eps2 = data(2:end,8);
delta_eps2(strcmp('n/a',data(2:end,8))) = {NaN};
delta_eps2 = cell2mat(delta_eps2);
delta_eps2 = delta_eps2(ind);
fr2 = data(2:end,9);
fr2(strcmp('n/a',data(2:end,9))) = {NaN};
fr2 = cell2mat(fr2);
fr2 = fr2(ind);
alpha2 = data(2:end,10);
alpha2(strcmp('n/a',data(2:end,10))) = {NaN};
alpha2 = cell2mat(alpha2);
alpha2 = alpha2(ind);
sigmaDC = data(2:end,11);
sigmaDC(strcmp('n/a',data(2:end,11))) = {0};
sigmaDC = cell2mat(sigmaDC);
sigmaDC = sigmaDC(ind);
Ea1 = data(2:end,12);
Ea1(strcmp('n/a',data(2:end,12))) = {NaN};
Ea1 = cell2mat(Ea1);
Ea1 = Ea1(ind);
Ea2 = data(2:end,13);
Ea2(strcmp('n/a',data(2:end,13))) = {NaN};
Ea2 = cell2mat(Ea2);
Ea2 = Ea2(ind);

k = 8.61733e-5; % eV/K


figure
%% Attenuation Rate (Temperature)
T = (223:273)'; % K
f = [100 500 1000 1500 2000]*1e6;
colors = plasma(length(f));
core = {'Vostok (Accreted 1)','Vostok (Meteoric)'};

% Cole-Cole
for i = 1:length(f)
    for n = 1:length(core)
        subplot(2,3,n)
        ind = find(strcmp(core(n),core_name));
        Na = zeros(length(T),length(ind));

        for m = 1:length(ind)

            Tref = T0(ind(m))+273.15;
            tau0 = 1./(2*pi*[fr1(ind(m)) fr2(ind(m))]);
            Ea = [Ea1(ind(m)) Ea2(ind(m))];
            tau = tau0.*exp((-Ea/k).*(1/Tref-1./T));
            delta_eps = [delta_eps1(ind(m)) delta_eps2(ind(m))];
            alpha = [alpha1(ind(m)) alpha2(ind(m))];

            if isnan(tau0(2)) || isnan(Ea(2)) || isnan(delta_eps(2)) || isnan(alpha(2))
                tau0(2) = [];
                Ea(2) = [];
                delta_eps(2) = [];
                alpha(2) = [];
            end


            eps_ice = colecole(eps_inf(ind(m)),delta_eps,tau,alpha,f(i),sigmaDC(ind(m)));
            [~, Na(:,m)] = EMalpha(eps_ice,f(i));
        end
        Na_mean = mean(Na,2);
        Na_std = std(Na,[],2);
        Na_min = Na_mean-Na_std;
        Na_max = Na_mean+Na_std;

        plot(T-273.15,Na_mean*1e3,'Color',colors(i,:));
        hold on
        patch([T; flipud(T)]-273.15,[Na_max*1e3; ...
            flipud(Na_min*1e3)],...
            colors(i,:),...
            'edgecolor','none',...
            'facealpha',0.1)
        hold on
    end

    % Matzler
    subplot(2,3,3)
    eps_ice = ice_matzler(T,f(i));
    [~, Na] = EMalpha(eps_ice,f(i));
    plot(T-273.15,Na*1e3,'Color',colors(i,:))
    hold on
end

subplot(2,3,1)
axis tight
ax(1) = gca;
grid on
ax(1).XLabel.String = 'Temperature, $T$ ($^{\circ}$C)';
ax(1).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(1).Title.String = 'Vostok (Accreted 1)';
ax(1).XLim = [-50 0];
ax(1).YLim = ylim2;


subplot(2,3,2)
axis tight
grid on
ax(2) = gca;
ax(2).XLabel.String = 'Temperature, $T$ ($^{\circ}$C)';
ax(2).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(2).Title.String = 'Vostok (Meteoric)';
ax(2).XLim = [-50 0];
ax(2).YLim = ylim1;

subplot(2,3,3)
axis tight
grid on
ax(3) = gca;
ax(3).XLabel.String = 'Temperature, $T$ ($^{\circ}$C)';
ax(3).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(3).Title.String = 'M$\ddot{\textrm{a}}$tzler';
ax(3).XLim = [-50 0];
ax(3).YLim = ylim1;

leg = legend(cellstr(num2str((f/1e9)')));
leg.Title.String = '$f$ (GHz)';
leg.Location = 'NorthWest';
leg.ItemTokenSize = [15 15];

%% Attenuation Rate (Frequency)
clear Na
f = logspace(-1,3).'*1e6;
T = [-40 -20 -10 -5 -1]+273.15;
core = {'Vostok (Accreted 1)','Vostok (Meteoric)'};
for i = 1:length(T)

    % Cole-Cole
    for n = 1:length(core)
        subplot(2,3,n+3)
        ind = find(strcmp(core(n),core_name));
        for m = 1:length(ind)

            Tref = T0(ind(m))+273.15;
            tau0 = 1./(2*pi*[fr1(ind(m)) fr2(ind(m))]);
            Ea = [Ea1(ind(m)) Ea2(ind(m))];
            tau = tau0.*exp((-Ea/k).*(1/Tref-1./T(i)));
            delta_eps = [delta_eps1(ind(m)) delta_eps2(ind(m))];
            alpha = [alpha1(ind(m)) alpha2(ind(m))];

            if isnan(tau0(2)) || isnan(Ea(2)) || isnan(delta_eps(2)) || isnan(alpha(2))
                tau0(2) = [];
                Ea(2) = [];
                delta_eps(2) = [];
                alpha(2) = [];
            end


            eps_ice = colecole(eps_inf(ind(m)),delta_eps,tau,alpha,f,sigmaDC(ind(m)));
            [~, Na(:,m)] = EMalpha(eps_ice,f);
        end
        Na_mean = mean(Na,2);
        Na_std = std(Na,[],2);
        Na_min = Na_mean-Na_std;
        Na_max = Na_mean+Na_std;
        semilogx(f,Na_mean*1e3,'Color',colors2(i,:));
        hold on
        patch([f; flipud(f)],[Na_max*1e3; ...
            flipud(Na_min*1e3)],...
            colors2(i,:),...
            'edgecolor','none',...
            'facealpha',0.1)
        hold on
    end

    % Matzler
    subplot(2,3,6)
    eps_ice = ice_matzler(T(i),f);
    [~, Na] = EMalpha(eps_ice,f);
    semilogx(f,Na*1e3,'Color',colors2(i,:))
    hold on
end

subplot(2,3,4)
ax(4) = gca;
grid on
ax(4).XMinorGrid = 'off';
axis tight
ax(4).XTick = 10.^(5:9);
ax(4).XLabel.String = 'Frequency, $f$';
ax(4).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(4).Title.String = 'Vostok (Accreted 1)';
ax(4).YLim = ylim2;

subplot(2,3,5)
ax(5) = gca;
axis tight
grid on
ax(5).XMinorGrid = 'off';
ax(5).XTick = 10.^(5:9);
ax(5).XLabel.String = 'Frequency, $f$';
ax(5).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(5).Title.String = 'Vostok (Meteoric)';
ax(5).YLim = ylim1;

subplot(2,3,6)
ax(6) = gca;
axis tight
grid on
ax(6).XMinorGrid = 'off';
ax(6).XTick = 10.^(5:9);
ax(6).XLabel.String = 'Frequency, $f$';
ax(6).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(6).Title.String = 'M$\ddot{\textrm{a}}$tzler';
ax(6).YLim = ylim1;

leg = legend(cellstr(num2str((T-273.15)')));
leg.Title.String = '$T$ ($^{\circ}$C)';
leg.Location = 'NorthWest';
leg.ItemTokenSize = [15 15];

%% Figure Formatting
%  [left bottom width height]

% scale
scale = 2;

% aspect ratio
AR = 1;

margin = 0.75;
h0 = 1*scale;
w0 = 1*AR*scale;
for n = 1:length(ax)
        % units
        ax(n).Units = 'inches';
        % height
        ax(n).Position(4) = h0;
        % width
        ax(n).Position(3) = w0;
end

% left
ax(1).Position(1) = margin;
ax(2).Position(1) = ax(1).Position(1)+ax(1).Position(3)+margin;
ax(3).Position(1) = ax(2).Position(1)+ax(2).Position(3)+margin;
ax(4).Position(1) = ax(1).Position(1);
ax(5).Position(1) = ax(2).Position(1);
ax(6).Position(1) = ax(3).Position(1);

% bottom
ax(4).Position(2) = 0.5*margin;
ax(5).Position(2) = ax(4).Position(2);
ax(6).Position(2) = ax(4).Position(2);
for n = 1:3
    ax(n).Position(2) = ax(4).Position(2)+ax(4).Position(4)+0.75*margin;
end


% labels
letters = {'a','b','c','d','e','f'};
for n = 1:length(letters)
        t(n) = annotation('textbox');
        t(n).String = letters{n};
        t(n).FontSize = fontsize+2;
        t(n).FontWeight = 'bold';
        t(n).Interpreter = 'tex';
        t(n).Units = 'inches';
        t(n).Position = [ax(n).Position(1)-0.75*margin ax(n).Position(2)+h0 0.1 0.1];
        t(n).EdgeColor = 'w';
        t(n).HorizontalAlignment = 'left';
        t(n).VerticalAlignment = 'bottom';
end

% Non-Zero Cole-Cole Distribution Parameter
t(1) = annotation('textbox');
t(1).String = 'Non-Zero Cole-Cole Distribution Parameter';
t(1).FontSize = fontsize+2;
t(1).FontWeight = 'bold';
t(1).Interpreter = 'tex';
t(1).Units = 'inches';
t(1).Position = [ax(1).Position(1) ax(1).Position(2)+ax(1).Position(4)+0.5*margin ax(2).Position(1)-ax(1).Position(1)+ax(2).Position(3) 0.1];
t(1).EdgeColor = 'w';
t(1).HorizontalAlignment = 'center';
t(1).VerticalAlignment = 'middle';

% IR Tail
t(2) = annotation('textbox');
t(2).String = 'Infrared Absorption';
t(2).FontSize = fontsize+2;
t(2).FontWeight = 'bold';
t(2).Interpreter = 'tex';
t(2).Units = 'inches';
t(2).Position = [ax(3).Position(1) ax(3).Position(2)+ax(3).Position(4)+0.5*margin ax(3).Position(3) 0.1];
t(2).EdgeColor = 'w';
t(2).HorizontalAlignment = 'center';
t(2).VerticalAlignment = 'middle';

f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 3*(margin+w0)+0.5*margin;
height = 2*(margin+h0);
f.Position(3:4) = [width height];
