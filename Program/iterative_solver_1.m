function [Scalar_flux, First_angular_moment, existing_angular_flux_positive, existing_angular_flux_negative, psi_positive, psi_negative] = iterative_solver_1(psi_boundary, L_plus_T_positive, L_plus_T_negative, Spatial_Material_Parameter, Angle_positive, Weight_positive, Angle_negative, Weight_negative, tau_z, Spatial_discretization_option)

% both directions are defined on the left boundary

error_z_0 = 1; % L2 error initilized. This is not compared with the true solution or finest mesh solution, but with different iterations. 
error_z_1 = 1;

[phi_z_0, phi_z_1] = Initial_gaess_phi_z_i(psi_boundary, L_plus_T_positive, L_plus_T_negative, Spatial_Material_Parameter, Angle_positive, Weight_positive, Angle_negative, Weight_negative, tau_z, Spatial_discretization_option); % Based on no scattering. z represents the index for the material length = length(Spatial_Material_Parameter.Thickness), i represent ith order of angular moment for the zeroth spatial moment

phi_z_0_positive = zeros(length(Spatial_Material_Parameter.Thickness), length(Angle_positive));
phi_z_1_positive = phi_z_0_positive;
phi_z_0_negative = phi_z_0_positive;
phi_z_1_negative = phi_z_0_positive;
psi_positive = phi_z_0_positive;
psi_negative = phi_z_0_positive;

existing_angular_flux_positive = zeros(1, length(Angle_positive));
existing_angular_flux_negative = existing_angular_flux_positive;

while error_z_0 > 10^(-6) && error_z_1 > 10^(-6)
    
    for i = 1: length(Angle_positive)
        %% positive directions
        [phi_z_0_positive(:, i), phi_z_1_positive(:, i), existing_angular_flux_positive(i), psi_positive(:, i)] = sovle_matrix_known_phi(phi_z_0, phi_z_1, Angle_positive(i), Spatial_Material_Parameter, tau_z(:, i), psi_boundary(i,1), L_plus_T_positive(:, :, i), 1, Weight_positive(i), Spatial_discretization_option);
        % for the existing_angular_flux_negative, it is actually not existing
        % but incoming. So, it might not be considered when presenting the data
        %% negative directions
        [phi_z_0_negative(:, i), phi_z_1_negative(:, i), existing_angular_flux_negative(i), psi_negative(:, i)] = sovle_matrix_known_phi(phi_z_0, phi_z_1, Angle_negative(i), Spatial_Material_Parameter,-tau_z(:, i), psi_boundary(i,2), L_plus_T_negative(:, :, i), 1, Weight_negative(i), Spatial_discretization_option);
    end
    
    phi_z_0_new = sum(phi_z_0_positive, 2) + sum(phi_z_0_negative, 2);
    phi_z_1_new = sum(phi_z_1_positive, 2) + sum(phi_z_1_negative, 2);
    
    error_z_0 = sqrt( sum( (phi_z_0_new - phi_z_0).^2 ) );
    error_z_1 = sqrt( sum( (phi_z_1_new - phi_z_1).^2 ) );
    
    phi_z_0 = phi_z_0_new;
    phi_z_1 = phi_z_1_new;
    
end

Scalar_flux = phi_z_0; % zeroth spatial moment and zeroth angular moment
First_angular_moment = phi_z_1; % zeroth spatial moment and first angular moment


end