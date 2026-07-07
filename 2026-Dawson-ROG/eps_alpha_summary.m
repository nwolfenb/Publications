% eps_alpha_summary.m
% Summary figure illustrating the dielectric response of ice
clear all; close all; clc

%% Path
addpath('..\..\IcyRF')

%% Parameters
eps_0 = 8.854e-12;
f = logspace(0,12,100).';
T = 250; % K
sigma_s = 1e-9;
eps_ice = ice_permittivity(T,f,sigma_s);
[alpha, Na] = EMalpha(eps_ice,f);
sigma = -imag(eps_ice).*(2*pi*eps_0*f);
sigma_hf = min(sigma(f>1e6));


figure
loglog(f,real(eps_ice),'r')
hold on
loglog(f,-imag(eps_ice),'b')
loglog(f,Na,'k')
loglog(f,sigma,'color',[0, 200, 0]/255)
loglog([1e6 max(f)],sigma_hf*ones(1,2),'--','color',[0, 200, 0]/255)

ax = gca;
ax.YLim = [1e-9 1e9];
ax.YTick = 10.^(-9:3:9);
ax.TickLabelInterpreter = 'tex';
ax.FontSize = 12;

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 7 5];
fig.Color = 'w';

