function [Macroscopic_Cross_Sections, Surface_Sources_info] = CrossSection(CrossSection_Source)

if strcmp( CrossSection_Source, 'HW')
%% Macroscopic Cross Sections

Material = {'Fuel';'Reflector';'Scatterer';'Absorber';'Isotropic';'Anisotropic';'Pure Absorber'};
Sigma_t = [1;2;2;10;0.1;0.1;10];
mu = [4;0;0;0;0;0;0];
Sigma_m = [0.2;0;0;0;0;0;0];
Sigma_s_0 = [0.6;1.15;1.27;2;0;0;0];
Sigma_s_1 = [0.1;-1;0;2;0;0;0];
Q_0 = [0;0;0;0;1;1;0];
Q_1 = [0;0;0;0;0;1;0];

Macroscopic_Cross_Sections = table(Sigma_t,mu, Sigma_m,Sigma_s_0,Sigma_s_1,Q_0,Q_1,'RowNames',Material);

%% Surface Sources

Source = {'Isotropic_surface';'Linear';'Discrete'};
Function = {'@(Omega_z)1'; '@(Omega_z)1+Omega_z'; '@(Omega_z)kroneckerDelta(Omega_z, 1)'}; % str2func(char(Surface_Sources{'Isotropic','Function'}))
Range_Omega_z = {'>0'; 'None'; 'None'};

Surface_Sources_info = table(Range_Omega_z, Function, 'RowNames', Source);
else
    warning('No CrossSection Data');
end
end