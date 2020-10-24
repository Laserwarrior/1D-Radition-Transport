function [phi_z_0, phi_z_1] = Initial_gaess_phi_z_i(psi_boundary, L_plus_T_positive, L_plus_T_negative, Spatial_Material_Parameter, Angle_positive, Weight_positive, Angle_negative, Weight_negative, tau_z, Spatial_discretization_option)

%% this is the gauss generator for angular moments of zeroth spatial moment at one angle
% based on only the Q^0 and Q^1 without scattering or multiplier
phi_z_0_positive = zeros(length(Spatial_Material_Parameter.Thickness), length(Angle_positive));
phi_z_1_positive = phi_z_0_positive;
phi_z_0_negative = phi_z_0_positive;
phi_z_1_negative = phi_z_0_positive;

for i = 1: length(Angle_positive)
    [phi_z_0_positive(:, i), phi_z_1_positive(:, i)] = sovle_matrix_not_known_phi(Angle_positive(i), Spatial_Material_Parameter, tau_z(:, i), psi_boundary(i,1), L_plus_T_positive(:, :, i), 1, Weight_positive(i), Spatial_discretization_option);
    [phi_z_0_negative(:, i), phi_z_1_negative(:, i)] = sovle_matrix_not_known_phi(Angle_negative(i), Spatial_Material_Parameter,-tau_z(:, i), psi_boundary(i,2), L_plus_T_negative(:, :, i), 1, Weight_negative(i), Spatial_discretization_option);  
end
    phi_z_0 = sum(phi_z_0_positive, 2) + sum(phi_z_0_negative, 2);
    phi_z_1 = sum(phi_z_1_positive, 2) + sum(phi_z_1_negative, 2);
    
end