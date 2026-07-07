% ColeCole_plot.m
% Cole-Cole plot for different models of relative permittivity
clear all; close all; clc

%% Path
addpath('..\..\IcyRF')

%% Defaults
fontsize = 8;
linewidth = 2;
colors = colororder;
markersize = 4;

%% Constants
f = logspace(0,12,1e3)'; % 1 Hz to 1 GHz
sigma = 0;
eps_inf = 3;
eps_s = [80 50];
delta_eps = eps_s-eps_inf;
fr1 = 5e3;
fr2 = 500;
fr = [fr1 fr2];
tau0 = 1./(2*pi*fr);
alpha1 = 0.1;
alpha2 = 0.05;
alpha = [alpha1 alpha2];

f0 = [1 1e3 10e3 100e3 1e6];
psym = {'o','s','^','d','p'};
frequencies = {'1 Hz','1 kHz','10 kHz','100 kHz','1 MHz'};

%% Models
figure
models = {'Cole-Cole (Multiple-Relaxations)','Cole-Cole (Single-Relaxation)','Debye (Multiple-Relaxations)','Debye (Single-Relaxation)'};
linestyle = {'-','-','-','-',':'};
for n = 1:length(models)
    if strcmp(models{n},'Cole-Cole (Multiple-Relaxations)')
        eps_ice  = colecole(eps_inf,delta_eps,tau0,alpha,f,sigma);

    elseif strcmp(models{n},'Cole-Cole (Single-Relaxation)')
        eps_ice  = colecole(eps_inf,delta_eps(1),tau0(1),alpha(1),f,sigma);

    elseif strcmp(models{n},'Debye (Multiple-Relaxations)')
        eps_ice  = colecole(eps_inf,delta_eps,tau0,[0 0],f,sigma);

    elseif strcmp(models{n},'Debye (Single-Relaxation)')
        eps_ice = debye(eps_s(1),eps_inf,tau0(1),f,sigma);

    else

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

    disp(models{n})
    eps_real = real(eps_ice);
    eps_imag = abs(imag(eps_ice));

    p(n) = plot(eps_real,eps_imag,'Color',colors(n,:),'LineWidth',linewidth,'LineStyle',linestyle{n});
    hold on

    for m = 1:length(f0)
        indf0 = find(f>=f0(m),1,'first');
        if n == length(models)
            p(m+n) = plot(eps_real(indf0),eps_imag(indf0),psym{m},...
                'Color','k','MarkerFaceColor','k','MarkerSize',markersize,'LineWidth',linewidth/2);
        else
            plot(eps_real(indf0),eps_imag(indf0),psym{m},...
                'MarkerFaceColor','k','Color','k','MarkerSize',markersize,'LineWidth',linewidth/2);
        end
    end

end

ax = gca;
grid on
ax.XLabel.String = 'Real Part, $\varepsilon^{\prime}$';
ax.YLabel.String = 'Imaginary Part, $\varepsilon^{\prime\prime}$';
leg = legend(p,[models,frequencies]);
leg.Location = 'EastOutside';
leg.ItemTokenSize = [15 15];
leg.FontSize = fontsize;


%% Annotations
x = (eps_s(1)+eps_inf)/2;
y = (eps_s(1)-eps_inf)/2;

% Vertical Line
plot(x*ones(1,2),[ax.YLim(1) ax.YLim(2)],'k:','LineWidth',linewidth/2)

% Horizontal Line
plot([ax.XLim(1) ax.XLim(2)],y*ones(1,2),'k:','LineWidth',linewidth/2)

% Labels
t1 = text(eps_inf+4, 1, '$\varepsilon_{\infty}$');
t1.HorizontalAlignment = 'left';
t1.VerticalAlignment = 'bottom';
t1.FontSize = fontsize;

t2 = text(eps_s(1)+4, 1, '$\varepsilon_{s}$');
t2.HorizontalAlignment = 'left';
t2.VerticalAlignment = 'bottom';
t2.FontSize = fontsize;

t4 = text(ax.XLim(2)+4,(eps_s(1)-eps_inf)/2, '$\displaystyle \frac{(\varepsilon_{s}-\varepsilon_{\infty})}{2}$');
t4.HorizontalAlignment = 'left';
t4.VerticalAlignment = 'bottom';
t4.FontSize = fontsize;

t5 = text((eps_s(1)+eps_inf)/2,ax.YLim(2)+1, '$\displaystyle \frac{(\varepsilon_{s}+\varepsilon_{\infty})}{2}$');
t5.HorizontalAlignment = 'center';
t5.VerticalAlignment = 'bottom';
t5.FontSize = fontsize;

%% Figure Formatting
%  [left bottom width height]

% scale
scale = 2;

% aspect ratio
AR = 1.5;

margin = 0.75;
h0 = 1*scale;
w0 = 1*AR*scale;
% units
ax.Units = 'inches';
% height
ax.Position(4) = h0;
% width
ax.Position(3) = w0;
% bottom
ax.Position(2) = 0.5*margin;
% left
ax.Position(1) = 0.75*margin;
% legend
leg.Units = 'inches';

f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = (margin+w0)+leg.Position(3)+margin/4;
height = (margin+h0);
f.Position(3:4) = [width height];

