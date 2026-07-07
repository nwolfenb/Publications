% att_low_loss.m
clear all; close all; clc

%% Path
addpath('..\..\IcyRF')

%% Defaults
colors = colororder;
fontsize = 8;
linewidth = 2;

%% Constants
f = 100e6;
eps_real = 3.17;
eps0 = 8.85418782e-12;
c = 3e8;

%% Loss Tangent

% Approximation
tand = logspace(-3,0,1e3);
tand_crit = 0.1;
term_full = (1+tand.^2).^(1/2);
term_approx = 1 + tand.^2/2;

% Conductivity to Loss Tangent
omega = 2*pi*f;
eps_imag = tand*eps_real;
sigma = eps_imag*omega*eps0;
eps_imag_crit = tand_crit*eps_real;
sigma_crit = eps_imag_crit*omega*eps0;

% Attenuation Rate
eps_r = eps_real+1j*eps_imag;
[alpha, Na] = EMalpha(eps_r,f);
att = Na*1e3;
att_approx = (10*log10(exp(1))/(3e8*eps0*sqrt(eps_real)))*sigma*1e3;

%% Figure
figure
% Approximation
subplot(1,3,1)
delta = term_approx-term_full;
loglog(tand,delta,'k','LineWidth',linewidth);
hold on
loglog([tand_crit tand_crit],...
    [min(delta) max(delta)],'k:','linewidth',1)
axis tight
grid on
ax(1) = gca;
ax(1).XMinorGrid = 'off';
ax(1).YMinorGrid = 'off';
ax(1).XTick = 10.^(-3:1);
ax(1).XLabel.String = 'Loss Tangent, $\tan{\delta}$';
ax(1).YLabel.String = '$(1+\tan^2\delta/2)-\sqrt{1+\tan^2\delta}$';
ax(1).Title.String = 'Taylor Expansion Error';

% Attenuation Rate
subplot(1,3,2)
delta = att_approx-att;
ind = find(delta<0);
loglog(tand,delta,'k','LineWidth',linewidth);
hold on
loglog(tand(ind),abs(delta(ind)),'k:','LineWidth',linewidth);
axis tight
loglog([tand_crit tand_crit],[min(abs(delta)) max(delta)],'k:','LineWidth',linewidth/2);

axis tight
grid on
ax(2) = gca;
ax(2).XMinorGrid = 'off';
ax(2).YMinorGrid = 'off';
ax(2).XTick = 10.^(-3:1);
ax(2).XLabel.String =  'Loss Tangent, $\tan{\delta}$';
ax(2).YLabel.String =  '$|\Delta N_{\alpha}|$ (dB/km)';
ax(2).Title.String = 'Attenuation Rate Error';
leg = legend('$\Delta N_{\alpha}>0$','$\Delta N_{\alpha}<0$');
leg.ItemTokenSize = [15 15];
leg.Location = 'NorthWest';

% Conductivity
subplot(1,3,3)
loglog(tand,sigma,'k','LineWidth',linewidth);
hold on
axis tight
loglog([tand_crit tand_crit],[min(sigma) max(sigma)],'k:','LineWidth',linewidth/2)

axis tight
grid on
ax(3) = gca;
ax(3).XMinorGrid = 'off';
ax(3).YMinorGrid = 'off';
ax(3).XTick = 10.^(-3:1);
ax(3).XLabel.String =  'Loss Tangent, $\tan{\delta}$';
ax(3).YLabel.String = 'Conductivity, $\sigma$ (S/m)';

%% Figure Formatting
%  [left bottom width height]

% scale
scale = 1.5;

% aspect ratio
AR = 1;

margin = 0.75;
h0 = 1*scale;
w0 = 1*AR*scale;
for p = 1:length(ax)
    % units
    ax(p).Units = 'inches';
    % height
    ax(p).Position(4) = h0;
    % width
    ax(p).Position(3) = w0;
    % bottom
    ax(p).Position(2) = 0.5*margin;
end

% left
ax(1).Position(1) = 0.75*margin;
for n = 2:length(ax)
    ax(n).Position(1) = ax(n-1).Position(1)+ax(n-1).Position(3)+margin;
end

% labels
letters = {'a','b','c'};
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


f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 3*(margin+w0);
height = (margin+h0);
f.Position(3:4) = [width height];