% Na_vs_T.m
clear all; close all; clc

%% Path
addpath('..\..\IcyRF')

%% Defaults
fontsize = 8;

%% Uncertainty
N = 1e3;

%% Constants
T = (223:273)'; % K
f0 = 100e6; % Hz
omega0 = 2*pi*f0;
ylim = [0 100];
eps0 = 8.854e-12;

%% Reference model
ref_model = 'Auty & Cole (1952)';

eps_inf_ref = 3.1;
eps_s_ref = static_permittivity(T,ref_model,N,false);
tau_ref = relaxation_time(T,ref_model,N,false);
eps_ice_ref = debye(eps_s_ref,eps_inf_ref,tau_ref,f0,0);
[~, Na_ref] = EMalpha(eps_ice_ref,f0);

%% Arrhenius
k = 8.617332e-5; % Boltzmann constant (eV/K)

% MacGregor et al. (2007), Table 1
Tr = 251; % Reference temperature (K)
Ea_ice = 0.55; % Pure ice (eV)
C_ice = 6.6; % Pure ice (uS/m)
Ea_ice_std = 0.05; % Pure ice (eV)
C_ice_std = 2.4; % Pure ice (uS/m)

N = 1e3;
Ea_ice_rand = Ea_ice + Ea_ice_std.*randn(1,N);
C_ice_rand = C_ice + C_ice_std.*randn(1,N);


subplot(2,2,2)
sigma = C_ice_rand.*exp(-(Ea_ice_rand/k).*((1./T)-(1/Tr)));
sigma_mat = mean(sigma,2) + [std(sigma,0,2) zeros(size(T)) -std(sigma,0,2)];
Na = 0.0009*sigma_mat*1e3;
p2(1) = plot(T-273.15,Na(:,2),'color','k');
hold on
patch([T; flipud(T)]-273.15,[Na(:,1); flipud(Na(:,3))],'k',...
    'edgecolor','none','facealpha',0.25)

% W97
Tr = -15+273.15; % C
Ea_ice = 0.58; % Pure ice (eV)
C_ice = 9; % Pure ice (uS/m)

sigma = C_ice.*exp(-(Ea_ice/k).*((1./T)-(1/Tr)));
Na = 0.0009*sigma*1e3;
p2(2) = plot(T-273.15,Na,'color',[240 71 15]/255);

% M07
Tr = -27+273.15; % Reference temperature (K)
Ea_ice = 0.51; % Pure ice (eV)
C_ice = 9.2; % Pure ice (uS/m)
Ea_ice_std = 0.01; % Pure ice (eV)
C_ice_std = 0.2; % Pure ice (uS/m)

N = 1e3;
Ea_ice_rand = Ea_ice + Ea_ice_std.*randn(1,N);
C_ice_rand = C_ice + C_ice_std.*randn(1,N);

sigma = C_ice_rand.*exp(-(Ea_ice_rand/k).*((1./T)-(1/Tr)));
sigma_mat = mean(sigma,2) + [std(sigma,0,2) zeros(size(T)) -std(sigma,0,2)];
Na = 0.0009*sigma_mat*1e3;
p2(3) = plot(T-273.15,Na(:,2),'color',[34 73 168]/255);
hold on
patch([T; flipud(T)]-273.15,[Na(:,1); flipud(Na(:,3))],[34 73 168]/255,'edgecolor','none','facealpha',0.25)

% % Check if this can ever exceed low loss approximation
% eps0 = 8.854e-12;
% omega = 2*pi*1e6;
% tand_max = max(1e-6*sigma_mat(:,1)/(omega*eps0*3.1));

% Johari & Charette
dataJC = readcell('Debye.xlsx','Sheet','Johari & Charette (1975)');
TJC = cell2mat(dataJC(2:end,3));
fJC = cell2mat(dataJC(2:end,5));
epsJC = cell2mat(dataJC(2:end,1));
tandJC = cell2mat(dataJC(2:end,2));
fJCu = unique(fJC);
psym = {'v','^'};
for n = 1:length(fJCu)
    ind = find(fJCu(n)==fJC & tandJC~=0);
    eps_real = epsJC(ind);
    tand = tandJC(ind);
    eps_imag = tand.*eps_real;
    eps_r = eps_real-1j*eps_imag;
    [~, Na] = EMalpha(eps_r,f0);
    p2(3+n) = plot(TJC(ind),Na*1e3,psym{n},'Color',[34 73 168]/255,'MarkerFaceColor',[34 73 168]/255);
end
% Reference
plot(T-273.15,Na_ref*1e3,'k--')


axis tight
ax(2) = gca;
ax(2).XLabel.String = 'Temperature, $T$ ($^{\circ}$C)';
ax(2).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(2).YLim = ylim;
ax(2).Layer = 'top';
leg(2) = legend(p2,'MacGregor et al. (2007)','W97','M07','35 MHz','60 MHz');
leg(2).ItemTokenSize = [15 15];
leg(2).Location = 'NorthWest';

grid on

%% Debye

% All models from Debye parameters plot
models1 = {'Auty & Cole (1952)',...
    'Worz & Cole (1969)',...
    'Gough & Davidson (1970)',...
    'Gough (1972)',...
    'Hasted (1973)',...
    'Johari & Charette (1975)',...
    'Johari & Jones (1975)',...
    'Johari & Jones (1976) - D2O',...
    'Kawada & Niinuma (1977)',...
    'Kawada (1978)',...
    'Kawada (1978) - parallel',...
    'Kawada (1978) - perpendicular',...
    'Johari & Jones (1978)',...
    'Johari & Jones (1978) - parallel',...
    'Johari & Jones (1978) - perpendicular',...
    'Johari & Whalley (1981)',...
    'Petrenko & Whitworth (1999)',...
    'Bittelli et al. (2004)',...
    'Sasaki et al. (2016) - Iha',...
    'Sasaki et al. (2016) - Ihb',...
    'Sasaki et al. (2016) - Ihc'};
colors1 = brewermap(length(models1),'spectral');

% Static Permittivity Models
eps_s_models = {'Auty & Cole (1952)',...
    'Worz & Cole (1969)',...
    'Gough & Davidson (1970)',...
    'Johari & Jones (1975)',...
    'Kawada & Niinuma (1977)',...
    'Kawada (1978)',...
    'Kawada (1978) - parallel',...
    'Kawada (1978) - perpendicular',...
    'Johari & Jones (1978)',...
    'Johari & Jones (1978) - parallel',...
    'Johari & Jones (1978) - perpendicular',...
    'Johari & Whalley (1981)'};

% Relaxation Time Models
tau_models = {'Auty & Cole (1952)',...
    'Kawada (1978)',...
    'Johari & Jones (1978)',...
    'Sasaki et al. (2016) - Iha',...
    'Sasaki et al. (2016) - Ihb',...
    'Sasaki et al. (2016) - Ihc'};
    % 'Bittelli et al. (2004)',... -> one of their published parameters is
    % clearly incorrect
    % 'Johari & Jones (1975)',... -> limited region of validity
    %  'Johari & Whalley (1981)',... -> super out-of-family
    %     'Johari & Jones (1976) - D2O',... -> D2O
    %     'Sasaki et al. (2016) - Ihb',... -> removing for clarity

    % HF Permittivity
    eps_inf = ice_gough(T-273.15);

    % Attenuation Rate
    % Uncertainty Bounds
    sigma = 0;
    for n = 1:length(tau_models)
        tau = relaxation_time(T,tau_models{n},N,false);
        if strcmp(tau_models{n},'Johari & Jones (1978)')
            % published uncertainty is so large it swamps plot
            tau = mean(tau,2);
        end

        Na = zeros(length(T),length(tau_models));
        for m = 1:length(eps_s_models)
            eps_s = static_permittivity(T,eps_s_models{m},false);
            eps_ice = debye(eps_s,eps_inf,tau,f0,sigma);

            [~, Na] = EMalpha(eps_ice,f0);
            Na = Na*1e3;
            Na_mean(:,m) = mean(Na,2);
            Na_std = std(Na,0,2);
            Na_min(:,m) = Na_mean(:,m)-Na_std;
            Na_max(:,m) = Na_mean(:,m)+Na_std;
        end
        Na_mean = mean(Na_mean,2);
        Na_min = min(Na_min,[],2);
        Na_max = max(Na_max,[],2);

        subplot(2,2,1)
        % plot(T-273.15,Na_mean,':','Color',colors1(strcmp(models1,tau_models{n}),:));
        patch([T; flipud(T)]-273.15,[Na_max; ...
            flipud(Na_min)],...
            colors1(strcmp(models1,tau_models{n}),:),...
            'edgecolor','none',...
            'facealpha',0.25)
        hold on
    end


% Mean Value
k = 0;
sigma = 0;
for m = 1:length(eps_s_models)
    eps_s = static_permittivity(T,eps_s_models{m},N,true);
    for n = 1:length(tau_models)
        tau = relaxation_time(T,tau_models{n},N,true);
        if strcmp(tau_models{n},'Johari & Jones (1978)')
            % published uncertainty is so large it swamps plot
            tau = mean(tau,2);
        end

        eps_ice = debye(eps_s,eps_inf,tau,f0,sigma);
        if ~all(isnan(eps_ice))
            k = k+1;

            [~, Na] = EMalpha(eps_ice,f0);

            Na = Na*1e3;
            Na_mean = mean(Na,2);
            Na_std = std(Na,0,2);
            Na_min = Na_mean-Na_std;
            Na_max = Na_mean+Na_std;

            p1(n) = plot(T-273.15,Na_mean,'Color',colors1(strcmp(models1,tau_models{n}),:));
            hold on
            ind = find(~isnan(Na_mean));
            % patch([T(ind); flipud(T(ind))]-273.15,[Na_max(ind); ...
            %     flipud(Na_min(ind))],...
            %     colors1(strcmp(models1,tau_models{n}),:),...
            %     'edgecolor','none',...
            %     'facealpha',0.25)
            % models{k,1} = [num2str(k),': ', eps_s_models{m},' + ',tau_models{n}];
            % Trange(k,:) = [min(T(~isnan(Na_mean))) max(T(~isnan(Na_mean)))];
            models{n} = tau_models{n};

        end
    end
end

% Reference
plot(T-273.15,Na_ref*1e3,'k--')

axis tight
ax(1) = gca;
ax(1).XLabel.String = 'Temperature, $T$ ($^{\circ}$C)';
ax(1).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(1).Title.String = 'Debye Model (Single-Relaxation)';
ax(1).YLim = ylim;
ax(1).Layer = 'top';
ind = cellfun(@(s) isnumeric(s) && isempty(s), models);
p1(ind) = [];
models(ind) = [];
models = cellfun(@(s) strrep(s, '&', '\&'), models, 'UniformOutput', false);


leg(1) = legend(p1,models);
leg(1).Location = 'NorthWest';
leg(1).ItemTokenSize = [15 15];
box on
grid on

%% Cole-Cole
clear Na

% Table 2, Stillman et al. (2013) - Acids
data_chem = readcell('Chemistry.xlsx');
core_name_chem = data_chem(2:end,1);
depth_chem = round(cell2mat(data_chem(2:end,2)));

% Table 2, Grimm et al. (2015)
data = readcell('Cole-Cole.xlsx');
core_name = data(2:end,1);
depth = round(cell2mat(data(2:end,2)));

% Isolate Shared Data
ind = find(ismember(core_name, core_name_chem));

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

core = unique(core_name);
colors = brewermap(length(core),'Set1');

k = 8.61733e-5; % eV/K
subplot(2,2,3)
% Min and Max Attenuation Rate
for n = 1:length(core)
    ind = find(strcmp(core(n),core_name));
    % Removing Samples from Anomalous Depths
    if strcmp(core(n),'Siple Dome') || strcmp(core(n),'GISP2') || strcmp(core(n),'Newall')
        ind(isnan(Ea2(ind))) = [];
    end
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


        eps_ice = colecole(eps_inf(ind(m)),delta_eps,tau,alpha,f0,sigmaDC(ind(m)));
        [~, Na(:,m)] = EMalpha(eps_ice,f0);

    end
    Na_mean = mean(Na,2);
    Na_std = std(Na,[],2);
    Na_min =  Na_mean-Na_std;
    Na_max = Na_mean+Na_std;

    patch([T; flipud(T)]-273.15,[Na_max; ...
        flipud(Na_min)]*1e3,...
        colors(n,:),...
        'edgecolor','none',...
        'facealpha',0.25)
    hold on
end

% Average Attenuation Rate
for n = 1:length(core)
    ind = find(strcmp(core(n),core_name));

    % Removing Samples from Anomalous Depths
    if strcmp(core(n),'Siple Dome') || strcmp(core(n),'GISP2') || strcmp(core(n),'Newall')
        ind(isnan(Ea2(ind))) = [];
    end
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


        eps_ice = colecole(eps_inf(ind(m)),delta_eps,tau,alpha,f0,sigmaDC(ind(m)));
        [~, Na(:,m)] = EMalpha(eps_ice,f0);

        % disp(['d = ',num2str(depth(ind(m))),' m'])
        % disp(['Na = ',num2str(min(Na(T>=round(Tref)))*1e3),' dB/km'])
    end
    Na_mean = mean(Na,2);
    p3(n) = plot(T-273.15,Na_mean*1e3,'Color',colors(n,:));
    hold on
end

% Reference
plot(T-273.15,Na_ref*1e3,'k--')

axis tight
ax(3) = gca;
ax(3).XLabel.String = 'Temperature, $T$ ($^{\circ}$C)';
ax(3).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(3).Title.String = 'Cole-Cole';
ax(3).YLim = ylim;
ax(3).Layer = 'top';
leg(3) = legend(p3,core);
leg(3).ItemTokenSize = [15 15];
leg(3).Location = 'NorthWest';
grid on
box on

%% Arrhenius w/ Impurities

% Chemistry
% Table 2, Grimm et al. (2015)
uMH = zeros(size(core));
uMCl = uMH;
uMNH4 = uMH;
for n = 1:length(core)
    ind = ismember(data_chem(:,1),core{n});
    uMH(n) = mean(cell2mat(data_chem(ind,3)));
    uMCl(n) = mean(cell2mat(data_chem(ind,4)));
    uMNH4(n) = mean(cell2mat(data_chem(ind,5)));
end


k = 8.617332e-5; % Boltzmann constant (eV/K)

% MacGregor et al. (2007), Table 1 & MacGregor et al. (2015), Table 2
Tr = 251; % Reference temperature (K)
Ea_ice = 0.55; % Pure ice (eV)
C_ice = 6.6; % Pure ice (uS/m)
Ea_Cl = 0.19; % Chloride (eV)
CM_Cl = 0.43; % Chloride (uS/m)
Ea_H = 0.20; % Acid (eV)
CM_H = 3.2; % Acid (uS/m)
Ea_NH4 = 0.23; % Ammonium (eV)
CM_NH4 = 0.8;  % Ammonium (uS/m)

Ea_ice_std = 0.05; % Pure ice (eV)
C_ice_std = 2.4; % Pure ice (uS/m)
Ea_Cl_std = 0.02; % Chloride (eV)
CM_Cl_std = 0.07; % Chloride (uS/m)
Ea_H_std = 0.04; % Acid (eV)
CM_H_std = 0.5; % Acid (uS/m)
Ea_NH4_std = 0;
CM_NH4_std = 0;

N = 1e3;
Ea_ice_rand = Ea_ice + Ea_ice_std.*randn(1,N);
C_ice_rand = C_ice + C_ice_std.*randn(1,N);
Ea_Cl_rand = Ea_Cl + Ea_Cl_std.*randn(1,N);
CM_Cl_rand = CM_Cl + CM_Cl_std.*randn(1,N);
Ea_H_rand = Ea_H + Ea_H_std.*randn(1,N);
CM_H_rand = CM_H + CM_H_std.*randn(1,N);
Ea_NH4_rand = Ea_NH4 + Ea_NH4_std.*randn(1,N);
CM_NH4_rand = CM_NH4 + CM_NH4_std.*randn(1,N);

subplot(2,2,4)

% cores = {'Vostok (Accreted)','Vostok (Meteoric)','Taylor','Siple Dome','GISP2'};
% % % Values from Eliza except for Vostok accreted acids
% uMCl = [1.72e-5 1.98e-6 7.33e-7 4.12e-6 1.15e-6];
% uMH = [3e-6 4.56e-7 1.57e-6 1.23e-6 7.44e-7];
% uMCl_std = [0 0 0 0 0];
% uMH_std = [4.4e-6 0 0 0 0];

% Uncertainty
for n = 1:length(core)

    C_Cl_rand = uMCl(n).*CM_Cl_rand;
    C_H_rand = uMH(n).*CM_H_rand;
    C_NH4_rand = uMNH4(n).*CM_NH4_rand;

    sigma = C_ice_rand.*exp(-(Ea_ice_rand/k).*((1./T)-(1/Tr)))+... % Pure
        C_Cl_rand.*exp(-(Ea_Cl_rand/k).*((1./T)-(1/Tr)))+... % Chloride
        C_H_rand.*exp(-(Ea_H_rand/k).*((1./T)-(1/Tr)))+... % Acid
        C_NH4_rand.*exp(-(Ea_NH4_rand/k).*((1./T)-(1/Tr)));
    sigma_mat = mean(sigma,2) + [std(sigma,0,2) zeros(size(T)) -std(sigma,0,2)];

    % Low loss approximation
    % Na = 0.0009*sigma_mat;

    % Full expression
    eps_p = ice_gough(T);
    eps_pp = 1e-6*sigma_mat./(omega0*eps0);
    eps_ice = eps_p - 1j*eps_pp;
    [~,Na] = EMalpha(eps_ice,f0);

    Na = 1e3*Na;
    patch([T; flipud(T)]-273.15,[Na(:,1); flipud(Na(:,3))],...
        colors(n,:),...
        'edgecolor','none','facealpha',0.25);
    hold on
end

% Mean
for n = 1:length(core)
    C_Cl_rand = uMCl(n)*CM_Cl_rand;
    C_H_rand = uMH(n)*CM_H_rand;

    sigma = C_ice_rand.*exp(-(Ea_ice_rand/k).*((1./T)-(1/Tr)))+... % Pure
        C_Cl_rand.*exp(-(Ea_Cl_rand/k).*((1./T)-(1/Tr)))+... % Chloride
        C_H_rand.*exp(-(Ea_H_rand/k).*((1./T)-(1/Tr))); % Acid
    sigma_mat = mean(sigma,2) + [std(sigma,0,2) zeros(size(T)) -std(sigma,0,2)];

    % low loss approximation
    % Na = 0.0009*sigma_mat;

    % full expression
    eps_p = ice_gough(T);
    eps_pp = 1e-6*sigma_mat./(omega0*eps0);
    eps_ice = eps_p - 1j*eps_pp;
    [~,Na] = EMalpha(eps_ice,f0);


    Na = 1e3*Na;
    p4(n) = plot(T-273.15,Na(:,2),'color',colors(n,:));
end

% Reference
plot(T-273.15,Na_ref*1e3,'k--')


axis tight
ax(4) = gca;
ax(4).XLabel.String = 'Temperature, $T$ ($^{\circ}$C)';
ax(4).YLabel.String = 'Attenuation Rate, $N_{\alpha}$ (dB/km)';
ax(4).YLim = ylim;
ax(4).Layer = 'top';
leg(4) = legend(p4,core);
leg(4).ItemTokenSize = [15 15];
leg(4).Location = 'NorthWest';
grid on
box on

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
ax(3).Position(1) = ax(1).Position(1);
ax(2).Position(1) = ax(1).Position(1)+ax(1).Position(3)+margin;
ax(4).Position(1) = ax(2).Position(1);

% bottom
ax(3).Position(2) = 0.5*margin;
ax(4).Position(2) = ax(3).Position(2);
for n = 1:2
    ax(n).Position(2) = ax(3).Position(2)+ax(3).Position(4)+0.75*margin;
end


% labels
letters = {'a','b','c','d'};
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

% Permittivity Model
t(1) = annotation('textbox');
t(1).String = 'Permittivity Model';
t(1).FontSize = fontsize+2;
t(1).FontWeight = 'bold';
t(1).Interpreter = 'tex';
t(1).Units = 'inches';
t(1).Position = [ax(1).Position(1) ax(1).Position(2)+ax(1).Position(4)+0.4*margin ax(1).Position(3) 0.1];
t(1).EdgeColor = 'w';
t(1).HorizontalAlignment = 'center';
t(1).VerticalAlignment = 'middle';

% Conductivity Model
t(2) = annotation('textbox');
t(2).String = 'HF Conductivity Model';
t(2).FontSize = fontsize+2;
t(2).FontWeight = 'bold';
t(2).Interpreter = 'tex';
t(2).Units = 'inches';
t(2).Position = [ax(2).Position(1) ax(2).Position(2)+ax(2).Position(4)+0.4*margin ax(2).Position(3) 0.1];
t(2).EdgeColor = 'w';
t(2).HorizontalAlignment = 'center';
t(2).VerticalAlignment = 'middle';

% Pure
t(3) = annotation('textbox');
t(3).String = 'Pure Laboratory Ice';
t(3).FontSize = fontsize+2;
t(3).FontWeight = 'bold';
t(3).Interpreter = 'tex';
t(3).Units = 'inches';
t(3).Position = [ax(1).Position(1)-0.7*margin ax(1).Position(2) w0 0.1];
t(3).EdgeColor = 'w';
t(3).HorizontalAlignment = 'center';
t(3).VerticalAlignment = 'middle';
t(3).Rotation = 90;


% Impure
t(4) = annotation('textbox');
t(4).String = 'Natural Ice';
t(4).FontSize = fontsize+2;
t(4).FontWeight = 'bold';
t(4).Interpreter = 'tex';
t(4).Units = 'inches';
t(4).Position = [ax(3).Position(1)-0.7*margin ax(3).Position(2) w0 0.1];
t(4).EdgeColor = 'w';
t(4).HorizontalAlignment = 'center';
t(4).VerticalAlignment = 'middle';
t(4).Rotation = 90;


f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 2*(margin+w0)+0.5*margin;
height = 2*(margin+h0);
f.Position(3:4) = [width height];
