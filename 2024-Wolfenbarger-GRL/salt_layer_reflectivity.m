% salt_layer_reflectivity.m 
% Used to generate the plot shown in Figure 2 in
%
% Wolfenbarger, N. S., Blankenship, D. D., Young, D. A., Scanlan, K. M.,
% Chivers, C. J., Findlay, D., Steinbruegge, G. B., Chan, K., Grima, C.,
% Soderlund, K. M., & Schroeder, D. M. (2024). Radar Characterization of
% Salt Layers in Europa’s Ice Shell as a Window into Critical Ice-Ocean
% Exchange Processes.

clear all; close all; clc

% Located at https://github.com/nwolfenb
addpath(genpath('..\..\IcyRF'))

%% Defaults
fontsize = 10;
linewidth = 2;

%% Radar Parameters
fc_vec = [9 9 60]*1e6; % center frequency (Hz)
BW_vec = [2.8 1 10]*1e6; % bandwidth (Hz)
fs_vec = [12 4.8 48]*1e6; % sampling frequency (Hz)
lw = [1 2 2];
ypos = [-10 -5 -15];
tau = 200e-6; % pulse length (s)

%% Constants
c = 3e8;
f = [9 60]*1e6; % Hz
lambda = c./f; % m

%% Constants
% MgSO4-11H2O
eps_salt = 4.9;
eps_ice = 3.1;
Vs_V_eut = 0.325;

%% Mixing
% Effective Permittivity
eps_e = eps_ice;
eps_i = eps_salt;
v = 2; % Bruggeman
eps_eff = mixing(eps_e,eps_i,Vs_V_eut,v);

[~, R_ideal] = EMcoef(eps_ice,eps_eff);

eps_r = [eps_ice eps_eff eps_ice];

%% Reflectivity
L = 1e3;

figure
% Extra subplots
subplot(2,2,1)
semilogx([7.9 7.9],[0 1],':','Color',[128,128,128]/255,'LineWidth',linewidth/2) % VHF
hold on
semilogx([79.1 79.1],[0 1],'k:','LineWidth',linewidth/2) % HF
xlim([1e-1 1e3])
ylim([0 1])
axis off
ax11 = gca;
ax11.XScale = 'log';


subplot(2,1,2)
p(1) = semilogx([1e-1 1e3],[R_ideal R_ideal],'k--','LineWidth',1);
hold on

h_min = zeros(size(fc_vec));
vres_thresh = h_min;
maxR = zeros(length(fc_vec),L);
for m = 1:length(fc_vec)
    lambda = c/fc_vec(m);
    fc = fc_vec(m);
    BW = BW_vec(m);
    fs = fs_vec(m);
    s = BW/tau; % chirp slope
    
    dt = 1/fs;
    N = floor(tau/dt);
    t = (-N/2:N/2-1)*dt;
    x = [exp(1i*(pi*s*t.^2 + 2*pi*fc*t)) zeros(1,1e3-length(t))];
    M = length(x);
    
    % FFT scaling note:
    % N = length of signal (not including zero-padding)
    % scale by 1/N to get units of dB (technically dB/sample?)
    % scale by dt to get units of dB/Hz
    % note that multiplying dt by df (to go from dB/Hz to dB)
    % dt * (1/(*dt)) = 1/N
    % also note that zero padding does not increase the energy of the
    % signal, which is why you scale by the non-zero padded length
    
    % frequency domain
    X = fft(x)/N;
    f = fs*(0:(M-1))/M;
    f(M/2+1:end) = f(M/2+1:end)-fs;
   
    % pulse-compression
    XX = X.*conj(X);
    
    %% Salt Layer (Thin Film)
    
    RF_dB = 20*log10(abs((sqrt(eps_r(1))-sqrt(eps_r(2)))./...
        (sqrt(eps_r(1))+sqrt(eps_r(2)))));
    
    n1 = sqrt(eps_r(1));
    n2 = sqrt(eps_r(2));
    n3 = sqrt(eps_r(3));
    
    r12 = (n2 - n1)/(n2+n1);
    r23 = (n3 - n2)/(n3+n2);
    
    while max(f)<fc % undersampled
        f = f+fs;
    end
    
    vres_thresh(m) = (c/n2)/(2*BW);
    t_crit = lambda/n2/4;
    
    l_vec = logspace(-1,3,L); % m
    h = l_vec;
    frac_lambda = h/lambda*n2;
    for n = 1:length(frac_lambda)
        beta = 2*pi*n2*h(n)./(c./f);
        R = (r12+r23*exp(2*1j*beta))./(1+r12*r23*exp(2*1j*beta));
        
        XXR = XX.*R;
        
        % note the scale factor applied to Mouginot et al. (2009)
        maxR(m,n) = 20*log10(max(abs(ifft(XXR)*N)));
    end
    [pks,locs] = findpeaks(maxR(m,:),h);
    h_min(m) = min(locs);

    if m == 1
        p(m+1) = semilogx(l_vec,maxR(m,:),'r','LineWidth',lw(m));
        hold on
        semilogx(vres_thresh(m)*ones(1,2),[-60 0],'r:','LineWidth',linewidth/2);
    elseif m == 2
        p(m+1) = semilogx(l_vec,maxR(m,:),'k','LineWidth',lw(m));
        hold on
        semilogx(vres_thresh(m)*ones(1,2),[-60 0],'k:','LineWidth',linewidth/2);
    else
        p(m+1) = semilogx(l_vec,maxR(m,:),'Color',[128,128,128]/255,'LineWidth',lw(m));
        hold on
        semilogx(vres_thresh(m)*ones(1,2),[-60 0],':','Color',[128,128,128]/255,'LineWidth',linewidth/2);
    end
end

axis tight
ax = gca;
ax.FontSize = fontsize;
ax.Children = [ax.Children(5); ax.Children(3); ax.Children(1); ax.Children(6); ax.Children(4); ax.Children(2); ax.Children(7)];
ax.XTick = [1e-1 1 1e1 1e2 1e3];
ax.XLabel.String = 'Salt Layer Thickness, $d_{\mathrm{salt}}$ (m)';
ax.XLabel.FontSize = fontsize;
ax.YLabel.String = 'Reflectivity, $R$ (dB)';
ax.YLabel.FontSize = fontsize;
ax.XLim = [1e-1 1e3];
ax.YLim = [-60 -20];
leg = legend(p,['$R=$ ',num2str(RF_dB,'%2.1f'),' dB'],'RIME','REASON HF','REASON VHF',...
'Location','SouthEast','FontSize',fontsize-2);
leg.ItemTokenSize(1) = leg.ItemTokenSize(2);

%% Below Range Resolution
text(7.9,-20,{'REASON     ';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','right','VerticalAlignment','Bottom')
text(7.9,-20,{'VHF';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')

text(28.2,-20,'RIME',...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')

text(79.1,-20,{'HF';''},...
    'FontSize',fontsize-2,'Interpreter','tex',...
    'HorizontalAlignment','center','VerticalAlignment','Bottom')


%% Figure Formatting
%  [left bottom width height]
M = 1;
N = 2;

% left
margin = 0.175;
ax.Position(1) = margin/N*N;
ax11.Position(1) = ax.Position(1);

% height
htot = 1-(margin/M+margin/2);
ax.Position(4) = htot/M;
ax11.Position(4) = margin/N/2;

% width
wtot = 1-(margin/N+3/4*margin);
ax.Position(3) = wtot;
ax11.Position(3) = ax.Position(3);

% bottom
ax.Position(2) = margin/M;
ax11.Position(2) = ax.Position(2)+ax.Position(4);


f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 3;
height = 3;
f.Position(3:4) = [width height];

