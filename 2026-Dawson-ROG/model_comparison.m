% model_comparison_v4.m
% Comparison of the relative permittivity, loss tangent, and attenuation
% rate for 
clear all; close all; clc

%% Path
addpath('..\..\IcyRF')

%% Defaults
fontsize = 8;
linewidth0 = 2;

%% Constants
f = logspace(0,12,1e3)'; % 1 Hz to 1 GHz
sigma0 = 1e-9;
eps_inf = 5;
eps_s = [75 100];
delta_eps = eps_s-eps_inf;
fr1 = 5e3;
fr2 = 500;
fr = [fr1 fr2];
tau0 = 1./(2*pi*fr);
alpha1 = 0.05;
alpha2 = 0.1;
alpha = [alpha1 alpha2];
eps0 = 8.854e-12;

%% Models
figure
models = {'Cole-Cole (Multiple-Relaxations)',...
    'Cole-Cole (Single-Relaxation)',...
    'Debye (Multiple-Relaxations)',...
    'Debye (Single-Relaxation)',...
    'M$\ddot{\textrm{a}}$tzler*'};
linewidth = linewidth0*[1 1 1/2 1/2 1];
linestyle = {'-','-','-','-',':'};
for n = 1:length(models)
    if strcmp(models{n},'Cole-Cole (Multiple-Relaxations)')
        eps_ice  = colecole(eps_inf,delta_eps,tau0,alpha,f,sigma0);

    elseif strcmp(models{n},'Cole-Cole (Single-Relaxation)')
        eps_ice  = colecole(eps_inf,delta_eps(1),tau0(1),alpha(1),f,sigma0);

    elseif strcmp(models{n},'Debye (Multiple-Relaxations)')
        eps_ice  = colecole(eps_inf,delta_eps,tau0,[0 0],f,sigma0);

    elseif strcmp(models{n},'Debye (Single-Relaxation)')
        eps_ice = debye(eps_s(1),eps_inf,tau0(1),f,sigma0);
    else

        eps_ice = debye(eps_s(1),eps_inf,tau0(1),f,sigma0);
        eps_real = real(eps_ice);
        eps_imag = -imag(eps_ice);

        alpha = eps_imag.*f/1e9;

        T = -10+273.15;

        f = f/1e9; % Hz to GHz
        B1 = 0.0207; % K/GHz
        b = 335; % K
        B2 = 1.16e-11; % GHz^-3
        dbeta = exp(-9.963+0.0372*(T-273.16));
        beta = (B1./T).*exp(b./T)./(exp(b./T)-1).^2+B2*f.^2+dbeta; % GHz^-1
        eps_imag = alpha./f+beta.*f;
        f = f*1e9; % GHz to Hz

        eps_ice = eps_real - 1j*eps_imag;
    end

    % Real Permittivity
    eps_real = real(eps_ice);
    subplot(5,1,1)
    p(n) =semilogx(f,eps_real,'LineWidth',linewidth(n),'LineStyle',linestyle{n});
    hold on

    % Imaginary Permittivity
    eps_imag = -imag(eps_ice);
    subplot(5,1,2)
    semilogx(f,eps_imag,'LineWidth',linewidth(n),'LineStyle',linestyle{n})
    hold on

    % Loss Tangent
    subplot(5,1,3)
    tand = eps_imag./eps_real;
    semilogx(f,tand,'LineWidth',linewidth(n),'LineStyle',linestyle{n})
    hold on

    % HF Conductivity
    subplot(5,1,4)
    omega = 2*pi*f;
    eps0 = 8.854e-12;
    sigma = omega*eps0.*eps_imag;
    loglog(f,sigma,'LineWidth',linewidth(n),'LineStyle',linestyle{n})
    hold on

    if strcmp(models{n},'Debye (Single-Relaxation)')
        sigma_inf = max(sigma);
    end

    % Attenuation Rate
    subplot(5,1,5)
    [~, Na] = EMalpha(eps_ice,f);
    semilogx(f,Na*1e3,'LineWidth',linewidth(n),'LineStyle',linestyle{n})
    hold on

end

% Real Permittivity
subplot(5,1,1)
ax(1) = gca;
grid on
ax(1).XLim = [min(f) max(f)];
ax(1).XTick = 10.^(0:3:12);
ax(1).YLim = [0 200];
ax(1).XTickLabel = {'$10^0$','$10^3$','$10^6$','$10^9$','$10^{12}$'};
ax(1).XLabel.String = 'Frequency, $f$ (Hz)';
ax(1).YLabel.String = '$\varepsilon^{\prime}$';

plot([fr1 fr1],ax(1).YLim,'k:','Linewidth',linewidth0/2)
plot([fr2 fr2],ax(1).YLim,'k:','Linewidth',linewidth0/2)

text(min(f),eps_s(1),'$\varepsilon_{s,1}$ ','HorizontalAlignment','Right');
text(min(f),eps_s(1)+eps_s(2),'$\varepsilon_{s,1}+\varepsilon_{s,2}$ ','HorizontalAlignment','Right');
text(max(f),eps_inf,'\ $\varepsilon_{\infty}$');


leg = legend(p,models);
leg.Location = 'NorthEast';
leg.ItemTokenSize = [15 15];

% Imaginary Permittivity
subplot(5,1,2)
ax(2) = gca;
grid on
ax(2).XLim = [min(f) max(f)];
ax(2).XTick = 10.^(0:3:12);
ax(2).XTickLabel = {'$10^0$','$10^3$','$10^6$','$10^9$','$10^{12}$'};
ax(2).XLabel.String = 'Frequency, $f$ (Hz)';
ax(2).YLabel.String = '$\varepsilon^{\prime\prime}$';
ax(2).YLim = [0 50];

plot([fr1 fr1],ax(2).YLim,'k:','Linewidth',linewidth0/2)
plot([fr2 fr2],ax(2).YLim,'k:','Linewidth',linewidth0/2)


% Loss Tangent
subplot(5,1,3)
ax(3) = gca;
grid on
ax(3).XLim = [min(f) max(f)];
ax(3).XTick = 10.^(0:3:12);
ax(3).YLim = [0 3];
ax(3).XTickLabel = {'$10^0$','$10^3$','$10^6$','$10^9$','$10^{12}$'};
ax(3).XLabel.String = 'Frequency, $f$ (Hz)';
ax(3).YLabel.String = '$\tan\delta$';


plot([fr1 fr1],ax(3).YLim,'k:','Linewidth',linewidth0/2)
plot([fr2 fr2],ax(3).YLim,'k:','Linewidth',linewidth0/2)

% Conductivity
subplot(5,1,4)
ax(4) = gca;
grid on
ax(4).XLim = [min(f) max(f)];
ax(4).YLim = [1e-10 1e-2];
ax(4).XTick = 10.^(0:3:12);
ax(4).XTickLabel = {'$10^0$','$10^3$','$10^6$','$10^9$','$10^{12}$'};
ax(4).XLabel.String = 'Frequency, $f$ (Hz)';
ax(4).YLabel.String = '$\sigma$ (S/m)';
ax(4).YMinorGrid = 'off';


loglog([fr1 fr1],ax(4).YLim,'k:','Linewidth',linewidth0/2)
loglog([fr2 fr2],ax(4).YLim,'k:','Linewidth',linewidth0/2)

text(min(f),sigma0,'$\sigma_s$ ','HorizontalAlignment','Right');
text(max(f),sigma_inf,'\ $\sigma_{\infty}$ \ ');

% Attenuation Rate
subplot(5,1,5)
ax(5) = gca;
grid on
ax(5).XLim = [min(f) max(f)];
ax(5).YLim = [0 100];
ax(5).XTick = 10.^(0:3:12);
ax(5).XTickLabel = {'$10^0$','$10^3$','$10^6$','$10^9$','$10^{12}$'};
ax(5).XLabel.String = 'Frequency, $f$ (Hz)';
ax(5).YLabel.String = '$N_{\alpha}$ (dB/km)';

plot([fr1 fr1],ax(5).YLim,'k:','Linewidth',linewidth0/2)
plot([fr2 fr2],ax(5).YLim,'k:','Linewidth',linewidth0/2)


%% Figure Formatting
%  [left bottom width height]

% scale
scale = 1.5;

% aspect ratio
AR = 1.75;

margin = 0.5;
h0 = 1*scale;
w0 = 1*AR*scale;
for p = 1:length(ax)
    % units
    ax(p).Units = 'inches';
    % left
    ax(p).Position(1) = margin;
    % height
    ax(p).Position(4) = h0;
    % width
    ax(p).Position(3) = w0;
end

% bottom
ax(5).Position(2) = 0.75*margin;
for n = fliplr(1:length(ax)-1)
    ax(n).Position(2) = ax(n+1).Position(2)+ax(n+1).Position(4)+0.5*margin;
end

% remove x labels
for n = 1:length(ax)-1
    ax(n).XLabel.String = [];
    ax(n).XTickLabel = [];
end


% labels
letters = {'a','b','c','d','e'};
for n = 1:length(letters)
    t(n) = annotation('textbox');
    t(n).String = letters{n};
    t(n).FontSize = fontsize+2;
    t(n).FontWeight = 'bold';
    t(n).Interpreter = 'tex';
    t(n).Units = 'inches';
    t(n).Position = [ax(n).Position(1)-margin ax(n).Position(2)+h0 0.1 0.1];
    t(n).EdgeColor = 'w';
    t(n).HorizontalAlignment = 'left';
    t(n).VerticalAlignment = 'middle';
end


f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = (margin+w0)+0.5*margin;
height = 5*(margin+h0)-1.75*margin;
f.Position(3:4) = [width height];


%% Insets

inset2 = axes;
copyobj(allchild(ax(2)),inset2);
inset2.Box = 'on';
inset2.Units = 'inches';
inset2.Position = ax(2).Position;
inset2.Position(1) = ax(2).Position(1)+0.5*w0;
inset2.Position(2) = ax(2).Position(2)+0.5*h0;
inset2.Position(3) = w0/2;
inset2.Position(4) = h0/2;
inset2.YTick = [];
inset2.XScale = 'log';
inset2.YScale = 'log';
inset2.XLim = [1e6 1e12];
% inset2.XTick = 10.^(6:3:9);
inset2.XTickLabel = [];
% grid(inset2,'on')
inset2.LineWidth = 0.5;
inset2.XColor = 0.5*ones(1,3);
inset2.YColor = 0.5*ones(1,3);

inset3 = axes;
copyobj(allchild(ax(3)),inset3);
inset3.Box = 'on';
inset3.Units = 'inches';
inset3.Position = ax(3).Position;
inset3.Position(1) = ax(3).Position(1)+0.5*w0;
inset3.Position(2) = ax(3).Position(2)+0.5*h0;
inset3.Position(3) = w0/2;
inset3.Position(4) = h0/2;
inset3.YTick = [];
inset3.XScale = 'log';
inset3.YScale = 'log';
inset3.XLim = [1e6 1e12];
% inset3.XTick = 10.^(6:3:9);
inset3.XTickLabel = [];
% grid(inset3,'on')
inset3.LineWidth = 0.5;
inset3.XColor = 0.5*ones(1,3);
inset3.YColor = 0.5*ones(1,3);



