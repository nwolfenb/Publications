% salt_layer_thickness_Sb0.m 
% Used to generate the plot shown in Figure 4 in
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

%% Defaults
fontsize = 10;
linewidth = 2;
interpreter = 'latex';
% brewermap available at MATLAB File Exchange
cmap = flipud(brewermap(256,'Spectral'));

%% Salt Layer
fn = '..\..\FreezingSimulations\PHREEQC\frezchem_ColdChem\MgSO4\MgSO4_1ppt.pqo';

[P,T,Sb,rho_b,sigma] = brine_properties_PHREEQC(fn);
Seut = max(Sb); % ppt
rhoeut = max(rho_b);

mm_MgSO4 = 120.37;
mm_MgSO411H2O = 318.535;
mv_MgSO411H2O = 207.44;

rhos = mm_MgSO411H2O/mv_MgSO411H2O; % g/cm^3

veut = 0.325;
S = linspace(0,Seut);
rhob = interp1(Sb,rho_b,S);
d0 = logspace(-1,3,1e3)';
d_salt = d0*S/1000.*(rhob/rhos)*mm_MgSO411H2O/mm_MgSO4*(1/veut);

%% Thickness

figure
subplot(2,2,3)

x = repmat(S,length(d0),1);
y = repmat(d0,1,length(S));

p = pcolor(x,y,d_salt);
p.LineStyle = 'none';
p.FaceColor = 'interp';
clim([1e-1 1e3])
axis tight

ax3 = gca;
ax3.ColorScale = 'log';
ax3.Colormap = cmap;
ax3.FontSize = fontsize;
ax3.TickLabelInterpreter = interpreter;
ax3.XLabel.String = 'Initial Brine Salinity, $S_{b,0}$ (ppt)';
ax3.XLabel.FontSize = fontsize;
ax3.XLabel.Interpreter = interpreter;
ax3.YLabel.String = 'Initial Reservoir Thickness, $d_0$ (m)';
ax3.YLabel.FontSize = fontsize;
ax3.YLabel.Interpreter = interpreter;
ax3.YScale = 'log';
ax3.YTick = [1e-1 1 1e1 1e2 1e3];
ax3.YTickLabel = {'$10^{-1}$','1','$10$','$10^{2}$','$10^{3}$'};
ax3.Layer = 'Top';

cb3 = colorbar;
cb3.FontSize = fontsize;
cb3.TickLabelInterpreter = 'latex';
cb3.Title.String = {'$d_{salt}$ (m)'};
cb3.Title.Interpreter = 'latex';
cb3.Limits = [1e-1 1e3];
cb3.Ticks = [1e-1 1 1e1 1e2 1e3];
cb3.TickLabels = {'$10^{-1}$','1','$10$','$10^{2}$','$10^{3}$'};

hold on
[C,h] = contour(x,y,d_salt,[1e-1 1 1e1 1e2 1e3],'LineColor','k');
clabel(C,h,'Interpreter','latex','FontSize',fontsize-2,'Color','k')
h.LabelSpacing = 1000;

%% Sb0
X = d_salt;
Y = repmat(d0,1,length(S));
V = repmat(S,length(d0),1);

x = X(:);
y = Y(:);
v = V(:);

xq = logspace(-1,3,1e3);
yq = logspace(-1,3,1e3);

Xq = repmat(xq,length(yq),1);
Yq = repmat(yq',1,length(xq));
F = scatteredInterpolant(x,y,v,'linear','none');
vq = F(Xq(:),Yq(:));
Vq = reshape(vq',length(yq),length(xq));

subplot(2,2,2)
semilogx([7.9 7.9],[0 1],'k:','LineWidth',linewidth/2) % VHF
hold on
semilogx([79.1 79.1],[0 1],'k:','LineWidth',linewidth/2) % HF
xlim([min(xq) max(xq)])
ylim([0 1])
axis off
ax2 = gca;
ax2.XScale = 'log';


subplot(2,2,4)
p = pcolor(Xq,Yq,Vq);
p.LineStyle = 'none';
p.FaceColor = 'interp';
clim([1 Seut])
axis tight

ax4 = gca;
ax4.ColorScale = 'log';
ax4.Colormap = cmap;
ax4.FontSize = fontsize;
ax4.TickLabelInterpreter = interpreter;
ax4.XLabel.String = 'Salt Layer Thickness, $d_{salt}$ (m)';
ax4.XLabel.FontSize = fontsize;
ax4.XLabel.Interpreter = interpreter;
ax4.YLabel.String = 'Initial Reservoir Thickness, $d_0$ (m)';
ax4.YLabel.FontSize = fontsize;
ax4.YLabel.Interpreter = interpreter;
ax4.XScale = 'log';
ax4.YScale = 'log';
ax4.XTick = [1e-1 1 1e1 1e2 1e3];
ax4.YTick = [1e-1 1 1e1 1e2 1e3];
ax4.XTickLabel = {'$10^{-1}$','1','$10$','$10^{2}$','$10^{3}$'};
ax4.XTickLabelRotation = 0;
ax4.YTickLabel = {'$10^{-1}$','1','$10$','$10^{2}$','$10^{3}$'};
ax4.YTickLabelRotation = 0;
ax4.Layer = 'Top';

cb4 = colorbar;
cb4.FontSize = fontsize;
cb4.TickLabelInterpreter = 'latex';
cb4.Title.String = {'$S_{b,0}$ (ppt)'};
cb4.Title.Interpreter = 'latex';
cb4.Limits = [1 Seut];
cb4.Ticks = [1 10 100];
cb4.TickLabels = {'1','10','100'};

keq = 0.07;
hold on
[C,h] = contour(Xq,Yq,Vq,round(Seut*keq)*ones(1,2),'LineColor','k');
clabel(C,h,'Interpreter','latex','FontSize',fontsize-2,'Color','k')
 h.LabelSpacing = 100;

% Below Range Resolution
patch([0.1 0.7 0.7 0.1],[0.1 0.1 1e3 1e3],0*ones(1,3),'EdgeColor','w')
text(exp((log(0.1)+log(0.7))/2),exp((log(0.1)+log(1e3))/2),'Undetectable (RIME, REASON)',...
    'FontSize',fontsize-1,'Rotation',90,'Interpreter','tex',...
    'HorizontalAlignment','Center','VerticalAlignment','Middle',...
    'FontWeight','Normal','Color','w')

patch([0.7 4.4 4.4 0.7],[0.1 0.1 1e3 1e3],0*ones(1,3),'EdgeColor','w')
text(exp((log(0.7)+log(4.4))/2),exp((log(0.1)+log(1e3))/2),'Undetectable (RIME, REASON HF)',...
    'FontSize',fontsize-1,'Rotation',90,'Interpreter','tex',...
    'HorizontalAlignment','Center','VerticalAlignment','Middle',...
    'FontWeight','Normal','Color','w')

plot([7.9 7.9],[0.1 1e3],'k:','LineWidth',linewidth/2)
text(7.9,1e3,{'REASON     ';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','right','VerticalAlignment','Bottom')
text(7.9,1e3,{'VHF';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')

plot([28.2 28.2],[0.1 1e3],'k:','LineWidth',linewidth/2)
text(28.2,1e3,'RIME',...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')

plot([79.1 79.1],[0.1 1e3],'k:','LineWidth',linewidth/2)
text(79.1,1e3,{'HF';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')

text(exp((log(7.9)+log(1e3))/2),exp((log(8.5)+log(1e3))/2),'Ocean Injection',...
    'FontSize',fontsize,'Rotation',45,'Interpreter','tex',...
    'HorizontalAlignment','Center','VerticalAlignment','Bottom',...
    'FontWeight','Bold','Color','w')


%% Figure Formatting
%  [left bottom width height]
M = 1;
N = 2;

% left
margin = 0.175;
ax3.Position(1) = margin/N;

% height
htot = 1-(margin/M+margin/2);
ax2.Position(4) = margin/N/2;
ax3.Position(4) = htot/M;
ax4.Position(4) = htot/M;


% width
wtot = 1-(margin/N+1.75*margin);
ax2.Position(3) = wtot/N;
ax3.Position(3) = wtot/N;
ax4.Position(3) = wtot/N;

% bottom
ax3.Position(2) = margin/M;
ax4.Position(2) = margin/M;
ax2.Position(2) = ax4.Position(2)+ax4.Position(4);

% left
ax4.Position(1) = ax3.Position(1)+ax3.Position(3)+2.25*margin/N; 
ax2.Position(1) = ax4.Position(1);

% labels
annotation('textbox',[ax3.Position(1)-margin/N ax3.Position(2)+ax3.Position(4)+margin/M/4/2 0.1 0.1],...
    'string','a','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')

annotation('textbox',[ax4.Position(1)-margin/N ax4.Position(2)+ax4.Position(4)+margin/M/4/2 0.1 0.1],...
    'string','b','fontsize',1.25*fontsize,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')

f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 7;
height = 3;
f.Position(3:4) = [width height];