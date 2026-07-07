% debye_parameters_uncertainty.m
% TO DO: Incorporate uncertanity into models, should be cake...hopefully

clear all; close all; clc

%% Path
addpath('..\..\IcyRF')

%% Defaults
fontsize = 8;

%% Enforce Validity
flag = 1;

%% Samples for Uncertainty Distribution
N = 1e3;

%% Temperature
T = (100:273)'; % K

%% Colormap (Debye Parameters)
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


%% Static Permittivity Models
% eps_s_models = {'Auty & Cole (1952)',...
%     'Worz & Cole (1969)',...
%     'Gough & Davidson (1970)',...
%     'Johari & Jones (1975)',...
%     'Kawada & Niinuma (1977)',...
%     'Kawada (1978)',...
%     'Johari & Jones (1978) - parallel',...
%     'Johari & Jones (1978) - perpendicular',...
%     'Johari & Jones (1978)',...
%     'Johari & Whalley (1981)'};

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
    'Johari & Whalley (1981)',...
    'Petrenko & Whitworth (1999)'};

%% Relaxation Time Models
% tau_models = {'Auty & Cole (1952)',...
%     'Johari & Jones (1975)',...
%     'Kawada (1978)',...
%     'Johari & Jones (1978) - parallel',...
%     'Johari & Jones (1978) - perpendicular',...
%     'Johari & Jones (1978)',...
%     'Bittelli et al. (2004)'};

tau_models = {'Auty & Cole (1952)',...
    'Worz & Cole (1969)',...
    'Gough & Davidson (1970)',...
    'Johari & Jones (1975)',...
    'Johari & Jones (1976) - D2O',...
    'Kawada (1978)',...
    'Johari & Jones (1978)',...
    'Johari & Whalley (1981)',...
    'Bittelli et al. (2004)',...
    'Sasaki et al. (2016) - Iha',...
    'Sasaki et al. (2016) - Ihb',...
    'Sasaki et al. (2016) - Ihc'};

%% HF Permittivity
source  = {'Gough (1972)','Auty & Cole (1952)',...
    'Hasted (1973)'};

figure
eps_inf = hf_permittivity(T,'Gough (1972)',N,flag);
eps_inf_mean = mean(eps_inf,2);
eps_inf_std = std(eps_inf,0,2);
eps_inf_min = eps_inf_mean-eps_inf_std;
eps_inf_max = eps_inf_mean+eps_inf_std;

subplot(2,3,2)
plot(T,eps_inf_mean,'Color',colors1(strcmp(models1,source{1}),:))
hold on
patch([T(~isnan(eps_inf_mean)); flipud(T(~isnan(eps_inf_mean)))],...
    [eps_inf_max(~isnan(eps_inf_mean));...
    flipud(eps_inf_min(~isnan(eps_inf_mean)))],...
    colors1(strcmp(models1,source{1}),:),...
    'edgecolor',colors1(strcmp(models1,source{1}),:),'facealpha',0.5)

for n = 1:length(models1)
    p(n) = plot([NaN NaN],[NaN NaN],'Color',colors1(n,:));
    hold on
end


for n = 2:length(source)
    % Hasted in Table 3 of Bittelli et al. (2004)
    tab = readtable('Debye.xlsx','Sheet',source{n});
    eps_inf = tab.eps_inf;
    T_inf = tab.T+273.15; % K
    plot(T_inf,eps_inf,'o','Color',colors1(strcmp(models1,source{n}),:),...
        'MarkerFaceColor',colors1(strcmp(models1,source{n}),:))
    hold on
end
axis tight
grid on
ax(2) = gca;
ax(2).XLabel.String = 'Temperature, $T$ (K)';
ax(2).YLabel.String = 'High Frequency Permittivity, $\varepsilon_{\infty}$';
ax(2).XLim = [100 273];
ax(2).XTick = 100:50:250;


labels = regexprep(models1, '&', '$\\&$');
labels = regexprep(labels, '- parallel', '$||$');
labels = regexprep(labels, '- perpendicular', '$\\perp$');
labels = regexprep(labels, 'D2O', 'D$_2$O');
leg = legend(p,labels);
leg.ItemTokenSize = [15 15];

%% Models
k = 0;
sigma = 0;
f = 100e6;
eps_inf = hf_permittivity(T,'Gough (1972)',N,flag);
for m = 1:length(eps_s_models)
    % Static Permittivity
    subplot(2,3,1)
    if strcmp(eps_s_models{m},'Auty & Cole (1952)')
        tab = readtable('Debye.xlsx','Sheet',eps_s_models{m});
        eps_s = tab.eps_s;
        T_s = tab.T+273.15; % K
        plot(T_s,eps_s,'o','Color',colors1(strcmp(models1,eps_s_models{m}),:),...
            'MarkerFaceColor',colors1(strcmp(models1,eps_s_models{m}),:))
        hold on

        eps_s = static_permittivity(T,eps_s_models{m},N,flag);
        plot(T,eps_s,'Color',colors1(strcmp(models1,eps_s_models{m}),:))

    elseif strcmp(eps_s_models{m},'Petrenko & Whitworth (1999)')
        % Table 3 of Bitelli et al. (2004)
        tab = readtable('Debye.xlsx','Sheet',eps_s_models{m});
        eps_s = tab.eps_s;
        T_s = tab.T+273.15; % K
        plot(T_s,eps_s,'o','Color',colors1(strcmp(models1,eps_s_models{m}),:),...
            'MarkerFaceColor',colors1(strcmp(models1,eps_s_models{m}),:))
        hold on
        continue
    else
        eps_s = static_permittivity(T,eps_s_models{m},N,flag);
        eps_s_mean = mean(eps_s,2);
        eps_s_std = std(eps_s,0,2);
        eps_s_min = eps_s_mean-eps_s_std;
        eps_s_max = eps_s_mean+eps_s_std;

        plot(T,eps_s_mean,'Color',colors1(strcmp(models1,eps_s_models{m}),:))
        hold on
        ind = find(~isnan(eps_s_mean));
        patch([T(ind); flipud(T(ind))],[eps_s_max(ind); ...
            flipud(eps_s_min(ind))],...
            colors1(strcmp(models1,eps_s_models{m}),:),...
            'edgecolor',colors1(strcmp(models1,eps_s_models{m}),:),...
            'facealpha',0.5)

    end

    % Relaxation Time
    for n = 1:length(tau_models)
        subplot(2,3,3)
        if strcmp(tau_models{n},'Auty & Cole (1952)')
            tau = relaxation_time(T,tau_models{n},N,flag);

            if m == 1
                semilogy(T,tau,'Color',colors1(strcmp(models1,tau_models{n}),:))
                hold on

                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                tau_tab = tab.tau;
                T_tab = tab.T+273.15; % K
                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
            end

        elseif strcmp(tau_models{n},'Worz & Cole (1969)')
            if m == 1
                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                tau_tab = tab.tau; % s
                T_tab = (1000./tab.x1000_T); % K

                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
                hold on
            end
            continue

        elseif strcmp(tau_models{n},'Gough & Davidson (1970)')
            if m == 1
                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                tau_tab = tab.tau; % s
                T_tab = (1000./tab.x1000_T); % K
                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
                hold on
            end
            continue

        elseif strcmp(tau_models{n},'Johari & Jones (1976) - D2O')
            tau = relaxation_time(T,tau_models{n},N,flag);

            if m == 1
                semilogy(T,tau,'Color',colors1(strcmp(models1,tau_models{n}),:))
                hold on

                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                tau_tab = tab.tauT./(1000./tab.x1000_T); % s
                T_tab = (1000./tab.x1000_T); % K
                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
            end

        elseif strcmp(tau_models{n},'Johari & Jones (1978)')
            tau = relaxation_time(T,tau_models{n},N,flag);
            if m == 1
                semilogy(T,mean(tau,2),'Color',colors1(strcmp(models1,tau_models{n}),:))
                hold on
            end

        elseif strcmp(tau_models{n},'Johari & Whalley (1981)')
            tau = relaxation_time(T,tau_models{n},N,flag);

            if m == 1
                semilogy(T,tau,'Color',colors1(strcmp(models1,tau_models{n}),:))
                hold on

                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                tau_tab = tab.tau; % s
                T_tab = (1000./tab.x1000_T); % K
                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
            end

        elseif strcmp(tau_models{n},'Bittelli et al. (2004)')
            if m == 1
                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                fr = tab.fr;
                w = 2*pi*fr;
                tau_tab = 1./w;
                T_tab = tab.T+273.15; % K
                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
            end
            continue
        elseif strcmp(tau_models{n},'Sasaki et al. (2016) - Iha')
            tau = relaxation_time(T,tau_models{n},N,flag);

            if m == 1
                semilogy(T,tau,'Color',colors1(strcmp(models1,tau_models{n}),:))
                hold on

                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                tau_tab = tab.tau; % s
                T_tab = (1000./tab.x1000_T); % K
                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
            end

        elseif strcmp(tau_models{n},'Sasaki et al. (2016) - Ihb')
            tau = relaxation_time(T,tau_models{n},N,flag);

            if m == 1
                semilogy(T,tau,'Color',colors1(strcmp(models1,tau_models{n}),:))
                hold on

                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                tau_tab = tab.tau; % s
                T_tab = (1000./tab.x1000_T); % K
                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
            end

        elseif strcmp(tau_models{n},'Sasaki et al. (2016) - Ihc')
            tau = relaxation_time(T,tau_models{n},N,flag);

            if m == 1
                semilogy(T,tau,'Color',colors1(strcmp(models1,tau_models{n}),:));
                hold on

                tab = readtable('Debye.xlsx','Sheet',tau_models{n});
                tau_tab = tab.tau; % s
                T_tab = (1000./tab.x1000_T); % K
                semilogy(T_tab,tau_tab,'o','Color',colors1(strcmp(models1,tau_models{n}),:),...
                    'MarkerFaceColor',colors1(strcmp(models1,tau_models{n}),:))
            end

        else
            tau = relaxation_time(T,tau_models{n},N,flag);

            if m == 1
                tau_mean = mean(tau,2);
                tau_std = std(tau,0,2);
                tau_min = tau_mean-tau_std;
                tau_min(tau_min<=0) = eps;
                tau_max = tau_mean+tau_std;

                semilogy(T,tau_mean,'Color',colors1(strcmp(models1,tau_models{n}),:))
                hold on
                ind = find(~isnan(tau_mean));
                patch([T(ind); flipud(T(ind))],[tau_max(ind); ...
                    flipud(tau_min(ind))],...
                    colors1(strcmp(models1,tau_models{n}),:),...
                    'edgecolor',colors1(strcmp(models1,tau_models{n}),:),...
                    'facealpha',0.5)
            end
        end

        % if size(eps_s,2)>1 && size(tau,2) == 1
        %     eps_ice = debye(eps_s(:,2),eps_inf,tau,f,sigma);
        % elseif size(eps_s,2) == 1 && size(tau,2)>1
        %     eps_ice = debye(eps_s,eps_inf,tau(:,2),f,sigma);
        % elseif size(eps_s,2)>1 && size(tau,2)>1
        %     eps_ice = debye(eps_s(:,2),eps_inf,tau(:,2),f,sigma);
        % else
        %     eps_ice = debye(eps_s,eps_inf,tau,f,sigma);
        % end
        eps_ice = debye(eps_s,eps_inf,tau,f,sigma);



        if ~all(isnan(eps_ice),'all')
            k = k+1;

            eps_real = real(eps_ice);
            eps_real_mean = mean(eps_real,2);
            eps_real_std = std(eps_real,0,2);
            eps_real_min = eps_real_mean-eps_real_std;
            eps_real_max = eps_real_mean+eps_real_std;


            eps_imag = -imag(eps_ice);
            eps_imag_mean = mean(eps_imag,2);
            eps_imag_std = std(eps_imag,0,2);
            eps_imag_min = eps_imag_mean-eps_imag_std;
            eps_imag_max = eps_imag_mean+eps_imag_std;


            subplot(2,3,5)
            if k == 1
                T_gough = linspace(2,273.15).';
                eps_real = hf_permittivity(T_gough,'Gough (1972)',N,flag);
                eps_real_mean = mean(eps_real,2);
                eps_real_std = std(eps_real,0,2);
                eps_real_min = eps_real_mean-eps_real_std;
                eps_real_max = eps_real_mean+eps_real_std;

                plot(T_gough,eps_real_mean,'k');
                hold on
                patch([T_gough(~isnan(eps_real_mean));...
                    flipud(T_gough(~isnan(eps_real_mean)))],...
                    [eps_real_max(~isnan(eps_real_mean)); ...
                    flipud(eps_real_min(~isnan(eps_real_mean)))],'k',...
                    'edgecolor','k','facealpha',0.5)
            end

            subplot(2,3,6)
            semilogy(T,eps_imag_mean,'Color',colors1(strcmp(models1,tau_models{n}),:))
            hold on

            models2{k,1} = [num2str(k),': ', eps_s_models{m},' + ',tau_models{n}];
            % Trange(k,:) = [min(T(~isnan(eps_ice_mean))) max(T(~isnan(eps_ice_mean)))];
        end
    end
end

subplot(2,3,1)
axis tight
grid on
ax(1) = gca;
ax(1).XLabel.String = 'Temperature, $T$ (K)';
ax(1).YLabel.String = 'Static Permittivity, $\varepsilon_s$';
ax(1).XLim = [100 273];
ax(1).XTick = 100:50:250;

subplot(2,3,3)
axis tight
grid on
ax(3) = gca;
ax(3).YMinorGrid = 'off';
ax(3).XLabel.String = 'Temperature, $T$ (K)';
ax(3).YLabel.String = 'Relaxation Time, $\tau$ (s)';
ax(3).YLim = [1e-5 1e3];
ax(3).YTick = [1e-5 1e-3 1e-1 1e1 1e3];
ax(3).XLim = [100 273];
ax(3).XTick = 100:50:250;

subplot(2,3,5)
axis tight
grid on
ax(5) = gca;
ax(5).XLabel.String = 'Temperature, $T$ (K)';
ax(5).YLabel.String = 'Real Part, $\varepsilon^{\prime}$';
ax(5).Title.String = ['$f=$ ',num2str(f/1e6),' MHz'];
ax(5).XLim = [100 273];
ax(5).XTick = 100:50:250;

subplot(2,3,6)
axis tight
grid on
ax(6) = gca;
ax(6).YMinorGrid = 'off';
ax(6).XLabel.String = 'Temperature, $T$ (K)';
ax(6).YLabel.String = 'Imaginary Part, $\varepsilon^{\prime\prime}$';
ax(6).Title.String = ['$f=$ ',num2str(f/1e6),' MHz'];
ax(6).YLim = [1e-10 1e-2];
ax(6).XLim = [100 273];
ax(6).XTick = 100:50:250;

%% Figure Formatting
%  [left bottom width height]

% scale
scale = 1.75;

% aspect ratio
AR = 1;

margin = 0.75;
h0 = 1*scale;
w0 = 1*AR*scale;
for n = 1:length(ax)
    if n ~=4
        % units
        ax(n).Units = 'inches';
        % height
        ax(n).Position(4) = h0;
        % width
        ax(n).Position(3) = w0;
    end
end

% left
ax(1).Position(1) = 0.75*margin;
for n = 2:3
    ax(n).Position(1) = ax(n-1).Position(1)+ax(n-1).Position(3)+margin;
    ax(n+3).Position(1) = ax(n).Position(1);
end

% bottom
ax(5).Position(2) = 1.5*margin;
ax(6).Position(2) = ax(5).Position(2);
for n = 1:3
    ax(n).Position(2) = ax(5).Position(2)+ax(5).Position(4)+0.75*margin;
end


% labels
letters = {'a','b','c',NaN,'d','e'};
for n = 1:length(letters)
    if ~isnan(letters{n})
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
end

%  [left bottom width height]
% legend
leg.Units = 'inches';
leg.NumColumns = 1;
leg.Position(1) = ax(1).Position(1);
leg.Position(2) = ax(5).Position(2) - (leg.Position(4)-ax(5).Position(4));


f = gcf;
f.Color = 'w';
f.Units = 'inches';
width = 3*(margin+w0);
height = 2*(margin+h0) + 0.75*margin;
f.Position(3:4) = [width height];





