% R_ice_salt.m 
% Used to generate the plot shown in Figure 1 in
%
% Wolfenbarger, N. S., Blankenship, D. D., Young, D. A., Scanlan, K. M.,
% Chivers, C. J., Findlay, D., Steinbruegge, G. B., Chan, K., Grima, C.,
% Soderlund, K. M., & Schroeder, D. M. (2024). Radar Characterization of
% Salt Layers in Europa’s Ice Shell as a Window into Critical Ice-Ocean
% Exchange Processes.

clear all; close all; clc

%% Add Paths
% Located at https://github.com/nwolfenb
addpath(genpath('..\..\FreezingSimulations'))
addpath(genpath('..\..\BrineVolumeFraction'))
addpath(genpath('..\..\IcyRF'))

%% Defaults
fontsize = 10;
linewidth = 2;
interpreter = 'latex';

%% Plot
figure
fdir = 'C:\Users\nswol\Dropbox\Research\Projects\PHREEQC\';
fn = 'frezchem_ColdChem\Binary\MgSO4_1ppt.pqo';

[~,T_liq,Sb_liq] = liquidus_PHREEQC([fdir,fn]);
Seut = max(Sb_liq);
S = linspace(0,Seut,1000);
T = linspace(50,160)-273.15;
xlims = [min(S) max(S)];

% Processes
name = {'Flash Freezing','Desalination','Cryoconcentration'};
Vpos = [0.25 0.75 0.75];
Halign = {'center','left','right'};
S_lim(1,:) = [0 155];
S_lim(2,:) = [0 12];
S_lim(3,:) = [155 167];
Hpos = [mean(S_lim(1,:)) min(S_lim(2,:)) max(S_lim(3,:))];


%% Labels
subplot(2,2,1)
xlim(xlims)
ylim([0 1])
axis off
ax11 = gca;

subplot(2,2,2)
for n = 1:length(S_lim)
    plot(S_lim(n,:),Vpos(n)*ones(1,2),'k','linewidth',4)
    text(Hpos(n),Vpos(n)+0.15,name{n},'fontsize',fontsize-2,...
        'fontweight','bold','VerticalAlignment','bottom','HorizontalAlignment',Halign{n},'interpreter','tex')
    hold on
end
hold on

xlim(xlims)
ylim([0 1])
axis off
ax22 = gca;

%% Salinity to Volume Fraction
[P, Tmat, Smat, Vi_V, Vb_V, Vs_V] = volume_fraction_PHREEQC(T,S,[fdir,fn]);
Vs_V_mean = mean(Vs_V,1)';
Vs_V_min = min(Vs_V,[],1)';
Vs_V_max = max(Vs_V,[],1)';

subplot(2,2,3)
patch([S'; flip(S')], [Vs_V_min; flip(Vs_V_max)], [1 1 1]*0.8, 'EdgeColor','none')
hold on
plot(S,Vs_V_mean,'k','LineWidth',linewidth);

ax1 = gca;
axis tight
ax1.FontSize = fontsize;
ax1.TickLabelInterpreter = interpreter;
ax1.XLabel.String = 'Bulk Ice Salinity, $S$ (ppt)';
ax1.XLabel.FontSize = fontsize;
ax1.XLabel.Interpreter = interpreter;
ax1.YLabel.String = 'Salt Hydrate Volume Fraction, $V_{ss}/V$';
ax1.YLabel.FontSize = fontsize;
ax1.YLabel.Interpreter = interpreter;
ax1.XMinorTick = 'on';
box on

%% Reflectivity
Rlim = -60; % dB (dynamic range)
[P, Tmat, Smat, Vi_V, Vb_V, Vs_V] = volume_fraction_PHREEQC(T,S,[fdir,fn]);
Vs_V = mean(Vs_V,1)';
v = 2; % Bruggeman

subplot(2,2,4)

% Minimum Bulk Ice Salinity Contrast
eps_ice = 3.0;
eps_salt = 5.2;

eps_eff = mixing(eps_ice,eps_salt,Vs_V,v);
[~, R_salt] = EMcoef(eps_eff,eps_eff');
R_salt(isinf(R_salt)) = NaN;


R_salt_vec = linspace(-60,max(R_salt,[],'all'),1000);
M = contourc(S,S,R_salt,R_salt_vec);

Mdata = contourdata(M); 
dS = zeros(1,length(Mdata)/2);
for n = 1:length(Mdata)/2
    xdata = [];
    ydata = [];
    for m = 1:2
        xdata = [xdata; Mdata(2*n-1+m-1).xdata];
        ydata = [ydata; Mdata(2*n-1+m-1).ydata];
    end
    dS(:,n)= min(abs(xdata-ydata));
end
Xmin = dS';
R_salt_vec_min = R_salt_vec';

% Maximum Bulk Ice Salinity Constrast
eps_ice = 3.2;
eps_salt = 4.6;

eps_eff = mixing(eps_ice,eps_salt,Vs_V,v);
[~, R_salt] = EMcoef(eps_eff,eps_eff');
R_salt(isinf(R_salt)) = NaN;


R_salt_vec = linspace(-60,max(R_salt,[],'all'),1000);
M = contourc(S,S,R_salt,R_salt_vec);

Mdata = contourdata(M); 
dS = zeros(1,length(Mdata)/2);
for n = 1:length(Mdata)/2
    xdata = [];
    ydata = [];
    for m = 1:2
        xdata = [xdata; Mdata(2*n-1+m-1).xdata];
        ydata = [ydata; Mdata(2*n-1+m-1).ydata];
    end
    dS(:,n)= max(abs(xdata-ydata));
end
Xmax = dS';
R_salt_vec_max = R_salt_vec';

patch([Xmin; flip(Xmax)], [R_salt_vec_min; flip(R_salt_vec_max)], [1 1 1]*0.9,'LineWidth',linewidth/4)
hold on


% Mean Bulk Ice Salinity Contrast
eps_ice = 3.1;
eps_salt = 4.9;

eps_eff = mixing(eps_ice,eps_salt,Vs_V,v);
[~, R_salt] = EMcoef(eps_eff,eps_eff');
R_salt(isinf(R_salt)) = NaN;


R_salt_vec = linspace(-60,max(R_salt,[],'all'),1000);
M = contourc(S,S,R_salt,R_salt_vec);

Mdata = contourdata(M); 
dS = zeros(1,length(Mdata)/2);
for n = 1:length(Mdata)/2
    xdata = [];
    ydata = [];
    for m = 1:2
        xdata = [xdata; Mdata(2*n-1+m-1).xdata];
        ydata = [ydata; Mdata(2*n-1+m-1).ydata];
    end
    dS(:,n)= mean(abs(xdata-ydata));
end
Xmean = dS';
R_salt_vec_mean = R_salt_vec';

plot(Xmean,R_salt_vec_mean,'k','LineWidth',linewidth)

axis tight
ax2 = gca;
ax2.FontSize = fontsize;
ax2.TickLabelInterpreter = interpreter;
ax2.XLabel.String = 'Bulk Ice Salinity Contrast, $\Delta S$ (ppt)';
ax2.XLabel.FontSize = fontsize;
ax2.XLabel.Interpreter = interpreter;
ax2.YLabel.String = 'Reflectivity, $R$ (dB)';
ax2.YLabel.FontSize = fontsize;
ax2.YLabel.Interpreter = interpreter;
ax2.YLim = [-60 -20];
ax2.XMinorTick = 'on';
box on



%% Minimum contrast that produces >-60 dB
M = contourc(S,S,R_salt,[-60 -60]);
Mdata = contourdata(M); 
dS = zeros(3,length(Mdata)/2);
for n = 1:length(Mdata)/2
    xdata = [];
    ydata = [];
    for m = 1:2
        xdata = [xdata; Mdata(2*n-1+m-1).xdata];
        ydata = [ydata; Mdata(2*n-1+m-1).ydata];
    end
    dS(:,n)= [min(abs(xdata-ydata)); mean(abs(xdata-ydata)); max(abs(xdata-ydata))];
end
disp(['min(dS) = ',num2str(dS(1,1)),' ppt'])

%% Figure Formatting
%  [left bottom width height]
M = 1;
N = 2;

% left
margin = 0.175;
ax1.Position(1) = margin/N;

% height
htot = 1-(margin/M+margin/2);
ax1.Position(4) = htot/M;
ax2.Position(4) = htot/M;
ax22.Position(4) = 1-htot/M-margin;

% width
wtot = 1-(margin/N+3/4*margin);
ax1.Position(3) = wtot/N;
ax2.Position(3) = wtot/N;
ax22.Position(3) = ax2.Position(3);

% left
ax2.Position(1) = ax1.Position(1)+ax1.Position(3)+1.25*margin/N; 
ax22.Position(1) =ax2.Position(1);

% bottom
ax2.Position(2) = margin/M-margin/4;
ax1.Position(2) = margin/M-margin/4;
ax22.Position(2) = ax2.Position(2) + ax2.Position(4);

% labels
annotation('textbox',[ax1.Position(1)-margin/N ax1.Position(2)+ax1.Position(4)+margin/M/4/2 0.1 0.1],...
    'string','a','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')

annotation('textbox',[ax2.Position(1)-margin/N ax2.Position(2)+ax2.Position(4)+margin/M/4/2 0.1 0.1],...
    'string','b','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')

f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 6;
height = 3;
f.Position(3:4) = [width height];
