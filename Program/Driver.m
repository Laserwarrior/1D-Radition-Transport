clear
clc

[CrossSection_DataSource, Surface_Sources_type, Material_distribution, Mesh_error_index, Mesh_error_num, Anlge_discretization_option, Spatial_discretization_option, Solving_method] = Input(); %% can be designed by manually input
[Macroscopic_Cross_Sections, Surface_Sources_info] = CrossSection(CrossSection_DataSource); 

if ~strcmp(Surface_Sources_type,'Vacuum')
    Surface_Sources = Surface_Sources_info(Surface_Sources_type,:);
    Surface_flag = 1;
else
    Surface_Sources = 'Vacuum';
    Surface_flag = 0;
end

[Angle_discretization_quadrature_set_name, Angle_discretization_quadrature_set_order] = Anlge_discretization_method(Anlge_discretization_option);

Existing_angular_flux_mesh_study = zeros(Mesh_error_num, str2double(Angle_discretization_quadrature_set_order)/2); % only consider even quadrature sets. Different column represents different directions. Different row coressponds to different mesh methods
% we cannot predefine the matrix to hold the information for scalar flux at
% different meshes, since the number of mesh is different 

%% analytic solution
% syms x
% psi_fun = @(x)exp(-10*(x)); % Pure absorber. may be modified to generate automatically based on the information in the future
% L2_Scalar_flux_Mesh_study_analytic = zeros(Mesh_error_num, 1);
%% solving the numerical solution

for i = 1:Mesh_error_num 
    
    Spatial_Material_Parameter = MeshGenerator(Material_distribution, Mesh_error_index, i, Macroscopic_Cross_Sections);
 
    [Scalar_flux, First_angular_moment, existing_angular_flux_positive, Angle_positive, psi_positive, Angle_negative, psi_negative, psi_boundary] = Solver(Surface_Sources, Surface_flag, Spatial_Material_Parameter, Anlge_discretization_option, Spatial_discretization_option, Solving_method);
    
    Existing_angular_flux_mesh_study(i,:) = existing_angular_flux_positive;
    
% %% used when analytic solution is known
%     Analytic_psi = psi_fun( Spatial_Material_Parameter.Material_right_boundary ); % based on the problem 2, here considers positive directions
%     Analytic_psi_index_minus_one = zeros(size(Analytic_psi));
%     Analytic_psi_index_minus_one(1) = psi_boundary(:,1)';
%     Analytic_psi_index_minus_one(2:end) = Analytic_psi(1:end-1); % based on the problem 2, here considers positive directions
%     phi_analytic = (Analytic_psi_index_minus_one - Analytic_psi)./(Spatial_Material_Parameter.Thickness.*Spatial_Material_Parameter.Sigma_t);
% %     discretized_analytic_error_plotter(Scalar_flux, phi_analytic, psi_positive, Analytic_psi, Spatial_Material_Parameter);
%     L2_Scalar_flux_Mesh_study_analytic(i) = sqrt( sum( abs( Scalar_flux - phi_analytic ).^2 ) ); % the scalar_flux should be a collum
    
end

% discretized_analytic_error_plotter(Scalar_flux, phi_analytic, psi_positive, Analytic_psi, Spatial_Material_Parameter);

Material_right_boundary = cumsum(Material_distribution.Thickness)*ones(1,2);
x_coordinate = Spatial_Material_Parameter.Material_right_boundary;
figure;
plot(x_coordinate, Scalar_flux,'LineStyle','-','Color','g','LineWidth',6); hold on
% plot(Material_right_boundary(1:end-1,:), [0, max(Scalar_flux)*1.2])
xlim([0, Material_right_boundary(end)])
ylim([0, max(Scalar_flux)*1.2])

figure;
plot(x_coordinate, First_angular_moment, 'LineStyle','-','Color','g','LineWidth',6); hold on
% plot(Material_right_boundary, [-max(First_angular_moment)*1.2, max(First_angular_moment)*1.2]);
xlim([0, Material_right_boundary(end)])
ylim([-max(First_angular_moment)*1.2, max(First_angular_moment)*1.2])

% L2_Scalar_flux_finest_Mesh_study = abs(Scalar_flux_mesh_study - Scalar_flux_mesh_study(end, :)); % cannot do it. mesh size change -> matrix change -> multigrid
L2_Existing_angular_flux_finest_Mesh_study = abs(Existing_angular_flux_mesh_study - Existing_angular_flux_mesh_study(end, :)); % Here for different angle, there is only one element. So, all the p-norm are all the same
Mesh_size_Mesh_study = Material_distribution.Mesh{Mesh_error_index};
h_Mesh_study = Material_distribution.Thickness(Mesh_error_index) ./ Mesh_size_Mesh_study;

% L2_Scalar_flux_finest_Mesh_except_last_one = L2_Scalar_flux_finest_Mesh_study(1:end-1);
L2_Existing_angular_flux_finest_Mesh_except_last_one = L2_Existing_angular_flux_finest_Mesh_study(1:end-1,1:end);
h_Mesh_study_except_last_one = h_Mesh_study(1:end-1);

% figure;plot(log10(h_Mesh_study_except_last_one), log10(L2_Scalar_flux_finest_Mesh_except_last_one)); % cannot do it. mesh size change -> matrix change -> multigrid
figure;plot(log10(h_Mesh_study_except_last_one), log10(L2_Existing_angular_flux_finest_Mesh_except_last_one));
% figure;plot(log10(h_Mesh_study),log10(L2_Scalar_flux_Mesh_study_analytic)); % use when analytic solution is known

log10_h_Mesh_study_except_last_one = log10(h_Mesh_study_except_last_one);
log10_L2_Existing_angular_flux_finest_Mesh_except_last_one = log10(L2_Existing_angular_flux_finest_Mesh_except_last_one);
% log10_L2_Scalar_flux_finest_Mesh_except_last_one = log10(L2_Scalar_flux_finest_Mesh_except_last_one); % cannot do it. mesh size change -> matrix change -> multigrid

%% used when analytic solution is known
% log10_h_Mesh_study = log10(h_Mesh_study);
% log10_L2_Scalar_flux_Mesh_study = log10(L2_Scalar_flux_Mesh_study_analytic);
