% Tb_data.m 
% Generates brightness temperature data for the publication
%
% Wolfenbarger, N. S., Broome, A. L., Schroeder, D. M., Ermakov, A. I.,
% Bolton, S. J., and Blankenship, D. D. (2025). Passive Microwave
% Radiometry and Active Radar Sounding as Complementary Tools for
% Geophysical Investigations of Icy Ocean Worlds.

clear all; close all; clc

%% Add Paths
% Located at https://github.com/nwolfenb
addpath(genpath('..\..\IcyRF'))
addpath(genpath('..\..\BrineVolumeFraction'))

%% Ideal
% Data for Figure 3 and Figures S5, S6 in Supporting Information S1
f = [0.6 1.2 2.5 5.0 10 22]'*1e9; % Hz
rho = 940; % kg/m^3
g = 1.315; % m/s^2
D = [1e-3 (1:1:100)]*1e3; % m
Ts = 100;
Tsky = 0;
rs = 0;
Tb = zeros(length(f),length(D));
Tb1 = Tb;
Tb2 = Tb;
Tb3 = Tb;
for m=1:length(f)
    for n = 1:length(D)
        P = rho*g*D(n);
        Tm = Tmelt(P);
        [rb, ~] = EMcoef(ice_permittivity(Tm-273.15,f(m),0),...
            water_permittivity(Tm-273.15,f(m)));

        z = (0:1:D(n))';
        T = Ts*(Tm/Ts).^(z./D(n));

        Z{m,n} = z;

        eps_ice = ice_permittivity(T-273.15,f(m),0);
        [Tb(m,n), Tb_z{m,n}, Tb1(m,n), Tb2(m,n), Tb3(m,n)] = brightness(T,z,eps_ice,rs,rb,f(m),Tsky);
    end
end

save('Tb.mat','Tb','f','D','Tb_z','Z','Tb1','Tb2','Tb3')

%% Log
% Data for Figure 2
f = [0.6 1.2 2.5 5.0 10 22]'*1e9; % Hz
rho = 940; % kg/m^3
g = 1.315; % m/s^2
D = logspace(0,5); % m
Ts = 100;
Tsky = 0;
rs = 0;
for m=1:length(f)
    for n = 1:length(D)
        P = rho*g*D(n);
        Tm = Tmelt(P);
        [rb, ~] = EMcoef(ice_permittivity(Tm-273.15,f(m),0),...
            water_permittivity(Tm-273.15,f(m)));

        z = (0:1:D(n))';
        T = Ts*(Tm/Ts).^(z./D(n));

        Z{m,n} = z;

        eps_ice = ice_permittivity(T-273.15,f(m),0);
        [Tb(m,n), Tb_z{m,n}, Tb1(m,n), Tb2(m,n), Tb3(m,n)] = brightness(T,z,eps_ice,rs,rb,f(m),Tsky);
    end
end

save('Tb_log.mat','Tb','f','D','Tb_z','Z','Tb1','Tb2','Tb3')

%% Convective Layer (260 K)
f = [0.6 1.2]'*1e9; % Hz
Tconv = 260;
fconv = 0:0.01:1;
D = [1e-3 (1:1:100)]*1e3; % m
Ts = 100;
Tsky = 0;
rs = 0;
Tb = zeros(length(f),length(D),length(fconv));
for q = 1:length(fconv)
    for m=1:length(f)
        for n = 1:length(D)

            % Base
            eps_ice = ice_permittivity(Tconv-273.15,f(m),0);

            [rb, ~] = EMcoef(eps_ice,...
                water_permittivity(Tconv-273.15,f(m)));

            z = (0:1:D(n))';
            dconv = D(n)*fconv(q);
            dcond = D(n) - dconv;
            T = zeros(size(z));
            if fconv(q) == 1
                T(z<=dcond) = Ts;
                T(z>dcond) = Tconv;
            else
                T(z<=dcond) = Ts*(Tconv/Ts).^(z(z<=dcond)./dcond);
                T(z>dcond) = Tconv;
            end

            % Ice Shell
            eps_ice = ice_permittivity(T-273.15,f(m),0);

            [Tb(m,n,q), ~, ~, ~, ~] = brightness(T,z,eps_ice,rs,rb,f(m),Tsky);
        end
    end
end

save('Tb_fconv_260.mat','Tb','f','D','fconv','Tconv')

%% Freezing Point Depression
f = [0.6 1.2]'*1e9; % Hz
deltaTm = 0:10:100;
rho = 940; % kg/m^3
g = 1.315; % m/s^2
D = [1e-3 (1:1:100)]*1e3; % m
Ts = 100;
Tsky = 0;
rs = 0;
Tb = zeros(length(f),length(D),length(deltaTm));
for q = 1:length(deltaTm)
    for m=1:length(f)
        for n = 1:length(D)
            P = rho*g*D(n);
            Tm = Tmelt(P)-deltaTm(q);

            % Base
            eps_ice = ice_permittivity(Tm-273.15,f(m),0);

            [rb, ~] = EMcoef(eps_ice,...
                water_permittivity(Tm-273.15,f(m)));

            z = (0:1:D(n))';
            T = Ts*(Tm/Ts).^(z./D(n));


            % Ice Shell
            eps_ice = ice_permittivity(T-273.15,f(m),0);

            [Tb(m,n,q), ~, ~, ~, ~] = brightness(T,z,eps_ice,rs,rb,f(m),Tsky);
        end
    end
end

save('Tb_deltaTm.mat','Tb','f','D','deltaTm')

%% Cl-Doped
f = [0.6 1.2]'*1e9; % Hz
uMCl = 0:30:300;
Tr = 251; % Reference temperature (K)
k = 8.617332e-5; % Boltzmann constant (eV/K)
CM_Cl = 0.5;
Ea_Cl = 0.19;
rho = 940; % kg/m^3
g = 1.315; % m/s^2
D = [1e-3 (1:1:100)]*1e3; % m
Ts = 100;
Tsky = 0;
rs = 0;
Tb = zeros(length(f),length(D),length(uMCl));
for q = 1:length(uMCl)
    C_Cl = uMCl(q)*CM_Cl;
    for m=1:length(f)
        for n = 1:length(D)
            P = rho*g*D(n);
            Tm = Tmelt(P);

            % Base
            eps_ice = ice_permittivity(Tm-273.15,f(m),0);

            [rb, ~] = EMcoef(eps_ice,...
                water_permittivity(Tm-273.15,f(m)));

            z = (0:1:D(n))';
            T = Ts*(Tm/Ts).^(z./D(n));
            
            % Ice Shell
            sigma_Cl = C_Cl.*exp(-(Ea_Cl/k).*((1./T)-(1/Tr)))/1e6;

            eps_ice = ice_permittivity(T-273.15,f(m),sigma_Cl);

            [Tb(m,n,q), ~, ~, ~, ~] = brightness(T,z,eps_ice,rs,rb,f(m),Tsky);
        end
    end
end

save('Tb_Cl.mat','Tb','f','D','uMCl')