% joint_constraints.m
% Generates Figure 6 for the publication 
% 
% Wolfenbarger, N. S., Broome, A. L., Schroeder, D. M., Ermakov, A. I.,
% Bolton, S. J., and Blankenship, D. D. (2025). Passive Microwave
% Radiometry and Active Radar Sounding as Complementary Tools for
% Geophysical Investigations of Icy Ocean Worlds.

clear all; close all; clc

%% Add Paths
% Located at https://github.com/nwolfenb
addpath(genpath('..\..\IcyRF'))

%% Defaults
fontsize = 8;
linewidth = 2;
interpreter = 'latex';

%% Parameters
d_obs = 15; % km
Na_obs = logspace(-1,1,3);
Ts = 100;

%% Freezing Point Depression
% Temperature Profile
p = 3;
subplot(3,2,p)
N = 1e4;
Tm = 250;
z = linspace(0,1,N)';
T = Ts*(Tm/Ts).^(z);

plot(T,z,'k','LineWidth',linewidth)
ax(p) = gca;
ax(p).YDir = 'reverse';
ax(p).FontSize = fontsize;
ax(p).TickLabelInterpreter = interpreter;
ax(p).XTick = [Ts Tm];
ax(p).XTickLabel = {'$100$';'$f(D)-\Delta T_f$'};
ax(p).XLabel.String = 'Temperature, $T$ (K)';
ax(p).XLabel.FontSize = fontsize;
ax(p).XLabel.Interpreter = interpreter;
ax(p).YLabel.String = 'Fractional Depth, $z/D$';
ax(p).YLabel.FontSize = fontsize;
ax(p).YLabel.Interpreter = interpreter;
ax(p).XLim = [Ts-50 270+50];
ax(p).YLim = [0 1];

t1 = text(Ts,0,'$T_s$');
t1.HorizontalAlignment='center';
t1.VerticalAlignment ='bottom';
t1.FontSize = fontsize;
t1.Interpreter = interpreter;

t2 = text(Tm,0,'$T_m$');
t2.HorizontalAlignment='center';
t2.VerticalAlignment ='bottom';
t2.FontSize = fontsize;
t2.Interpreter = interpreter;

% Joint Constraints (15 km)
p = 4;
subplot(3,2,p)
load('Tb_deltaTm.mat')
delta_Tb = squeeze(Tb(1,:,:)-Tb(2,:,:));
x = repmat(D'/1e3,1,length(deltaTm));
y = repmat(deltaTm,length(D),1);


[C,h] = contour(y.',x.',delta_Tb.',[10:5:40],'LineWidth',linewidth,'LineColor',[0 112 192]/255);
ylim([0 30])
clabel(C,h,'Interpreter',interpreter,'FontSize',fontsize,'Color',[0 112 192]/255,...
       'LabelSpacing',45)
hold on

f = 9e6;
D = 1e3;
dz = 1;
z = (0:dz:D)';
Tm0 = 273;
att = zeros(length(z),length(deltaTm));
for q = 1:length(deltaTm)
    Tm = Tm0-deltaTm(q);

    T = Ts*(Tm/Ts).^(z/D);
    eps_ice = ice_permittivity(T-273.15,f,0);

    [alpha, Na] = EMalpha(eps_ice,f);
    att(:,q) = 2e3*Na;

end

x = repmat(deltaTm,length(z),1);
y = repmat(z/1e3,1,length(deltaTm));

fill([0 100 100 0],[0 0 d_obs d_obs],'w','EdgeColor','none')

plot(deltaTm,d_obs*ones(size(deltaTm)),'r','LineWidth',linewidth/2);
hold on

[C,h] = contour(x,d_obs./y,att,Na_obs,'LineColor','r',...
    'LineWidth',linewidth);
clabel(C,h,'Interpreter',interpreter,'FontSize',fontsize,'Color','r',...
        'LabelSpacing',65);
text(0,d_obs,['$z_{echo}=$ ',num2str(d_obs),' km '],...
    'FontSize',fontsize,'Color','r',...
    'HorizontalAlignment','Left','VerticalAlignment','Top','Interpreter',interpreter)


ax(p) = gca;
ax(p).FontSize = fontsize;
ax(p).TickLabelInterpreter = interpreter;
ax(p).XLabel.String = 'Freezing Point Depression, $\Delta T_f$ (K)';
ax(p).XLabel.FontSize = fontsize;
ax(p).XLabel.Interpreter = interpreter;
ax(p).XLim = [0 100];
ax(p).XTick = 0:20:100;
ax(p).YLabel.String = 'Ice Shell Thickness, $D$ (km)';
ax(p).YLabel.FontSize = fontsize;
ax(p).YLabel.Interpreter = interpreter;
ax(p).YLim = [0 30];
ax(p).YTick = 0:10:30;
ax(p).Title.String = 'Joint Constraints';
ax(p).Title.FontSize = fontsize;
ax(p).Title.Interpreter = interpreter;

%% Convective Layer
% Temperature Profile
p = 1;
subplot(3,2,p)
Tm = 260;
N = 1e4/2;
z = linspace(0,1/2,N)';
T = Ts*(Tm/Ts).^(z/(1/2));
z = [z; linspace(1/2,1,N)'];
z = flipud(z);
T = [T; ones(size(T))*Tm];

plot(T,z,'k','LineWidth',linewidth)
ax(p) = gca;
ax(p).FontSize = fontsize;
ax(p).TickLabelInterpreter = interpreter;
ax(p).XTick = [Ts Tm];
ax(p).XTickLabel = {'$100$';'$260$'};
ax(p).XLabel.String = 'Temperature, $T$ (K)';
ax(p).XLabel.FontSize = fontsize;
ax(p).XLabel.Interpreter = interpreter;
ax(p).YLabel.String = 'Fractional Depth, $z/D$';
ax(p).YLabel.FontSize = fontsize;
ax(p).YLabel.Interpreter = interpreter;
ax(p).XLim = [Ts-50 270+50];
ax(p).YLim = [0 1];
ax(p).YTick = [0 1];
ax(p).YTickLabel = cellstr(num2str(([1 0])'));

t1 = text(Ts,1,'$T_s$');
t1.HorizontalAlignment='center';
t1.VerticalAlignment ='bottom';
t1.FontSize = fontsize;
t1.Interpreter = interpreter;

t2 = text(Tm,1,'$T_m$');
t2.HorizontalAlignment='center';
t2.VerticalAlignment ='bottom';
t2.FontSize = fontsize;
t2.Interpreter = interpreter;

delta_x = 15;
delta_y = 0.025;

ar = annotation('doublearrow');
ar.Parent = ax(p);
ar.X = [Tm-delta_x Tm-delta_x];
ar.Y = [0+delta_y 0.5-delta_y];
ar.Head1Length = 5;
ar.Head2Length = 5;
ar.Head1Width = 5;
ar.Head2Width = 5;
ar.Head1Style = 'plain';
ar.Head2Style = 'plain';
ar.LineWidth = linewidth/2;

txt = text(Tm-delta_x,0.25,'$d_{conv}/D$ ');
txt.FontSize = fontsize;
txt.HorizontalAlignment = 'right';
txt.Interpreter = interpreter;

% Joint Constraints (15 km)
p = 2;
subplot(3,2,p)
load('Tb_fconv_260.mat')
delta_Tb = squeeze(Tb(1,:,:)-Tb(2,:,:));
x = repmat(D'/1e3,1,length(fconv));
y = repmat(fconv,length(D),1);

[C,h] = contour(y.',x.',delta_Tb.',[10:10:40],'LineWidth',linewidth,'LineColor',[0 112 192]/255);
ylim([0 30])
clabel(C,h,'Interpreter',interpreter,'FontSize',fontsize,'Color',[0 112 192]/255,...
       'LabelSpacing',70)
hold on

f = 9e6;
D = 1e3;
dz = 1;
z = (0:dz:D)';
Ts = 100;
att = zeros(length(z),length(fconv));
for q = 1:length(fconv)
    
    T = zeros(size(z));
    dconv = D*fconv(q);
    dcond = D - dconv;
    T(z<=dcond) = Ts*(Tconv/Ts).^(z(z<=dcond)./dcond);
    T(z>dcond) = Tconv;

    eps_ice = ice_permittivity(T-273.15,f,0);

    [alpha, Na] = EMalpha(eps_ice,f);
    att(:,q) = 2e3*Na;

end

x = repmat(fconv,length(z),1);
y = repmat(z/1e3,1,length(fconv));


fill([0 1 1 0],[0 0 d_obs d_obs],'w','EdgeColor','none')

plot(fconv,d_obs*ones(size(fconv)),'r','LineWidth',linewidth/2);
hold on

[C,h] = contour(x,d_obs./y,att,Na_obs,'LineColor','r',...
    'LineWidth',linewidth);
ylim([0 30])
clabel(C,h,'Interpreter',interpreter,'FontSize',fontsize,'Color','r',...
        'LabelSpacing',70)

text(0,d_obs,['$z_{echo}=$ ',num2str(d_obs),' km '],...
    'FontSize',fontsize,'Color','r',...
    'HorizontalAlignment','Left','VerticalAlignment','Top','Interpreter',interpreter)


ax(p) = gca;
ax(p).FontSize = fontsize;
ax(p).TickLabelInterpreter = interpreter;
ax(p).XLabel.String = 'Convective Fraction, $d_{conv}/D$';
ax(p).XLabel.FontSize = fontsize;
ax(p).XLabel.Interpreter = interpreter;
ax(p).XLim = [0 1];
ax(p).XTick = 0:0.2:1;
ax(p).YLabel.String = 'Ice Shell Thickness, $D$ (km)';
ax(p).YLabel.FontSize = fontsize;
ax(p).YLabel.Interpreter = interpreter;
ax(p).YLim = [0 30];
ax(p).YTick = 0:10:30;
ax(p).Title.String = 'Joint Constraints';
ax(p).Title.FontSize = fontsize;
ax(p).Title.Interpreter = interpreter;

%% Chloride-Doped Ice
% Temperature Profile
p = 5;
subplot(3,2,p)
N = 1e4;
Tm = 270;
z = linspace(0,1,N)';
T = Ts*(Tm/Ts).^(z);

plot(T,z,'k','LineWidth',linewidth)
ax(p) = gca;
ax(p).YDir = 'reverse';
ax(p).FontSize = fontsize;
ax(p).TickLabelInterpreter = interpreter;
ax(p).XTick = [Ts Tm];
ax(p).XTickLabel = {'$100$';'$f(D)$'};
ax(p).XLabel.String = 'Temperature, $T$ (K)';
ax(p).XLabel.FontSize = fontsize;
ax(p).XLabel.Interpreter = interpreter;
ax(p).YLabel.String = 'Fractional Depth, $z/D$';
ax(p).YLabel.FontSize = fontsize;
ax(p).YLabel.Interpreter = interpreter;
ax(p).XLim = [Ts-50 270+50];
ax(p).YLim = [0 1];

t1 = text(Ts,0,'$T_s$');
t1.HorizontalAlignment='center';
t1.VerticalAlignment ='bottom';
t1.FontSize = fontsize;
t1.Interpreter = interpreter;

t2 = text(Tm,0,'$T_m$');
t2.HorizontalAlignment='center';
t2.VerticalAlignment ='bottom';
t2.FontSize = fontsize;
t2.Interpreter = interpreter;

% Joint Constraints (15 km)
p = 6;
subplot(3,2,p)
load('Tb_Cl.mat')
delta_Tb = squeeze(Tb(1,:,:)-Tb(2,:,:));
x = repmat(D'/1e3,1,length(uMCl));
y = repmat(uMCl,length(D),1);

[C,h] = contour(y.',x.',delta_Tb.',[10:5:40],'LineWidth',linewidth,'LineColor',[0 112 192]/255);
ylim([0 30])
clabel(C,h,'Interpreter',interpreter,'FontSize',fontsize,'Color',[0 112 192]/255,...
       'LabelSpacing',55)
hold on

f = 9e6;
Tr = 251;
k = 8.617332e-5;
CM_Cl = 0.5;
Ea_Cl = 0.19;
D = 1e3;
dz = 1;
z = (0:dz:D)';
Tm = 273;
att = zeros(length(z),length(uMCl));
for q = 1:length(uMCl)
    C_Cl = uMCl(q)*CM_Cl;

    T = Ts*(Tm/Ts).^(z./D);

    sigma_Cl = C_Cl.*exp(-(Ea_Cl/k).*((1./T)-(1/Tr)))/1e6;

    eps_ice = ice_permittivity(T-273.15,f,sigma_Cl);

    [alpha, Na] = EMalpha(eps_ice,f);
    att(:,q) = 2e3*Na;

end

x = repmat(uMCl,length(z),1);
y = repmat(z/1e3,1,length(uMCl));

fill([0 300 300 0],[0 0 d_obs d_obs ],'w','EdgeColor','none')
plot(uMCl,d_obs*ones(size(uMCl)),'r','LineWidth',linewidth/2);
hold on
[C,h] = contour(x,d_obs./y,att,Na_obs,'LineColor','r',...
    'LineWidth',linewidth);
ylim([0 30])
clabel(C,h,'Interpreter',interpreter,'FontSize',fontsize,'Color','r',...
        'LabelSpacing',50)
text(0,d_obs,['$z_{echo}=$ ',num2str(d_obs),' km '],...
    'FontSize',fontsize,'Color','r',...
    'HorizontalAlignment','Left','VerticalAlignment','Top','Interpreter',interpreter)

ax(p) = gca;
ax(p).FontSize = fontsize;
ax(p).TickLabelInterpreter = interpreter;
ax(p).XLabel.String = 'Ice Chlorinity, [Cl$^{-}$] ($\mu$M)';
ax(p).XLabel.FontSize = fontsize;
ax(p).XLabel.Interpreter = interpreter;
ax(p).XLim = [0 300];
ax(p).XTick = 0:50:300;
ax(p).YLabel.String = 'Ice Shell Thickness, $D$ (km)';
ax(p).YLabel.FontSize = fontsize;
ax(p).YLabel.Interpreter = interpreter;
ax(p).YLim = [0 30];
ax(p).YTick = 0:10:30;
ax(p).Title.String = 'Joint Constraints';
ax(p).Title.FontSize = fontsize;
ax(p).Title.Interpreter = interpreter;

%% Figure Formatting
%  [left bottom width height]

% units
for p = 1:length(ax)
ax(p).Units = 'centimeters';
end

% left
margin = 1;
ax(1).Position(1) = 1.5*margin;
ax(3).Position(1) = 1.5*margin;
ax(5).Position(1) = 1.5*margin;

% height
h0 = 4.25;
for p = 1:length(ax)
ax(p).Position(4) = h0;
end

% width
w0 = 4.25;
for p = 1:length(ax)
ax(p).Position(3) = w0;
end

% bottom
ax(5).Position(2) = margin;
ax(6).Position(2) = margin;
ax(3).Position(2) = ax(5).Position(2)+ax(5).Position(4)+1.5*margin;
ax(4).Position(2) = ax(6).Position(2)+ax(6).Position(4)+1.5*margin;
ax(1).Position(2) = ax(3).Position(2)+ax(3).Position(4)+1.5*margin;
ax(2).Position(2) = ax(4).Position(2)+ax(4).Position(4)+1.5*margin;

% left
ax(2).Position(1) = ax(1).Position(1)+ax(1).Position(3)+1.25*margin;
ax(4).Position(1) = ax(3).Position(1)+ax(3).Position(3)+1.25*margin;
ax(6).Position(1) = ax(5).Position(1)+ax(5).Position(3)+1.25*margin;

% labels
letters = {'a','b','c','d','e','f'};
for n = 1:length(letters)
    t(n) = annotation('textbox');
    t(n).String = letters{n};
    t(n).FontSize = fontsize+2;
    t(n).FontWeight = 'bold';
    t(n).Interpreter = 'tex';
    t(n).Units = 'centimeters';
    t(n).Position = [ax(n).Position(1)-margin ax(n).Position(2)+h0 0.1 0.1];
    t(n).EdgeColor = 'w';
    t(n).HorizontalAlignment = 'left';
    t(n).VerticalAlignment = 'bottom';
end

% Freezing Point Depression
n = 8;
t(n) = annotation('textbox');
t(n).String = 'Freezing Point Depression';
t(n).FontSize = fontsize+2;
t(n).FontWeight = 'bold';
t(n).Interpreter = 'tex';
t(n).Units = 'centimeters';
t(n).Position = [0 ax(3).Position(2)-0.5 h0+2*0.5 0.1];
t(n).EdgeColor = 'w';
t(n).HorizontalAlignment = 'center';
t(n).VerticalAlignment = 'top';
t(n).Rotation = 90;

% Convective Layer
n = 7;
t(n) = annotation('textbox');
t(n).String = 'Convective Layer';
t(n).FontSize = fontsize+2;
t(n).FontWeight = 'bold';
t(n).Interpreter = 'tex';
t(n).Units = 'centimeters';
t(n).Position = [0 ax(1).Position(2) h0 0.1];
t(n).EdgeColor = 'w';
t(n).HorizontalAlignment = 'center';
t(n).VerticalAlignment = 'top';
t(n).Rotation = 90;

% Chloride-Doped Ice
n = 9;
t(n) = annotation('textbox');
t(n).String = 'Chloride-Doped Ice';
t(n).FontSize = fontsize+2;
t(n).FontWeight = 'bold';
t(n).Interpreter = 'tex';
t(n).Units = 'centimeters';
t(n).Position = [0 ax(5).Position(2) h0 0.1];
t(n).EdgeColor = 'w';
t(n).HorizontalAlignment = 'center';
t(n).VerticalAlignment = 'top';
t(n).Rotation = 90;

f = gcf;
f.Color = 'w';
f.Units = 'centimeters';
width = 2*(1.25*margin+w0)+margin;
height = 3*(1.5*margin+h0);
f.Position(3:4) = [width height];