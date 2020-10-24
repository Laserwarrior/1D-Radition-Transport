function [CrossSection_DataSource, Surface_Sources_type, Material_distribution, Mesh_error_index, Mesh_error_num, Anlge_discretization_option, Spatial_discretization_option, Solving_method] = Input()

%% data used by the CrossSections. Here, we don't need to have the address, because we input the data manually
CrossSection_DataSource = 'HW';


%% Surface source type: boundary condition; The hw only has one side nontrival boundary condition
Surface_Sources_type = 'Isotropic_surface'; % left boundary condition. The right boundary condition is Vacuum for this problem. can be chosen from 'Discrete', 'Isotropic_surface', 'Linear', 'Vacuum'

%% spatial distribution of materials

%% Condition for Problem 4
% Material = {'Absorber';'Reflector'; 'Isotropic'; 'Scatterer'; 'Fuel'; 'Scatterer'; 'Isotropic'; 'Reflector'; 'Absorber'};
% Thickness = [   1;         1;         0.2;         0.1;        0.05;      0.1;         0.2;          1;         1]; % unit: cm
% Mesh = cell(length(Material), 1); % used for the error analysis. Now only code for one specific region
% 
% Mesh_error_index = 1; % index for which region to do the error analysis
% Mesh{Mesh_error_index} = [100];
% Mesh_error_num = length(Mesh{Mesh_error_index}); % number of different meshes to do the analysis
% 
% % Mesh{2} = [2]; % need to change according to the Mesh_error_index
% 
% for i = 2: length(Material)
%     Mesh{i} = Thickness(i)*100;
% end
% Material_distribution = table(Material, Thickness, Mesh);

%% Condition for pure absorber. Problem 2
Material = {'Scatterer'};
Thickness = [1]; % unit: cmD
Mesh = cell(length(Thickness), 1); % used for the error analysis. Now only code for one specific region

Mesh_error_index = 1; % index for which region to do the error analysis
% Mesh{Mesh_error_index} = [200, 400, 800, 1600];
Mesh{Mesh_error_index} = [200];
Mesh_error_num = length(Mesh{Mesh_error_index}); % number of different meshes to do the analysis

Material_distribution = table(Material, Thickness, Mesh);

%% Angle discretization options. Only use even angles 

% Anlge_discretization_option = 'Discrete:S2';

Anlge_discretization_option = 'LSGQ:S20'; % 2, 4, 6, 8, 12, 20

%% Spatial discretization

Spatial_discretization_option = 'Diamond Difference'; %'Step Characteristics', 'Diamond Difference'

%% Solving method

Solving_method = 'Transport_Sweep';


end
