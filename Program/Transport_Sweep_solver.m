function [Scalar_flux, First_angular_moment, existing_angular_flux_positive, Angle_positive, psi_positive, Angle_negative, psi_negative, psi_boundary] = Transport_Sweep_solver(Surface_Sources, Surface_flag, Spatial_Material_Parameter, Anlge_discretization_option, Spatial_discretization_option)

[Angle_discretization_quadrature_set_name, Angle_discretization_quadrature_set_order] = Anlge_discretization_method(Anlge_discretization_option);
[Angle, Weight] =  Angle_Weight_generator(Angle_discretization_quadrature_set_name, Angle_discretization_quadrature_set_order); % Angle will always be (- +), (- +)

Angle_positive = Angle(2:2:end); % number of Angle needs to be even
Weight_positive = Weight(2:2:end);
Angle_negative = Angle(1:2:end-1);
Weight_negative = Weight(1:2:end-1);

[psi_boundary, boundary_flag] = psi_boundary_fun(Surface_flag, Surface_Sources, Angle_positive, Angle_negative, Angle);

% boundary_flag : 1, Discrete boundary condition and Linear boundary condition. All the directions are
% defined on the left boundary. Because it is dirac, only omega_z = 1 would
% contribution. For both negative and positive omega_z,
% sweep from left to right. However, for the HW 1, since we only have pure
% absorber, it is not possible to build up from the material. Only need to
% sweep from the left to right.
% psi_boundary = [phi_0_positive, phi_0_negative]

% boundary_flag : 2, Isotropic and vacuum condition. Positive directions
% are defined on the left boundary, and negative directions are defined on
% the right boundary. For positive directions, sweep from left to right.
% For the negative directions, sweep from right to left.
% psi_boundary = [phi_0_positive, phi_N_negative]

% The solving process depends on the boundary flag
% For the existing angular flux, only positive directions need to be considered

tau_z =  Spatial_Material_Parameter.Sigma_t .* Spatial_Material_Parameter.Thickness * (1./Angle_positive)'; % each collumn represents different directions. each row represents different positions. This is absolute value 
Diagonal_ones =  eye(length(Spatial_Material_Parameter.Thickness));

L_plus_T_positive = zeros(length(Spatial_Material_Parameter.Thickness), length(Spatial_Material_Parameter.Thickness), length(Angle_positive));
L_plus_T_negative = zeros(length(Spatial_Material_Parameter.Thickness), length(Spatial_Material_Parameter.Thickness), length(Angle_positive));
        
switch boundary_flag
    case 1
        switch Spatial_discretization_option
            case 'Step Characteristics'
                for j = 1: length(Angle_positive) % use different pages to store the information for different angles
                    L_plus_T_positive(:, :, j) = Diagonal_ones;
                    L_plus_T_negative(:, :, j) = Diagonal_ones;
                    for i = 2: length(Spatial_Material_Parameter.Thickness)
                        L_plus_T_positive(i ,i-1, j) = - exp(-tau_z(i, j)); 
                        L_plus_T_negative(i ,i-1, j) = - exp( tau_z(i, j)); 
                    end
                end
            case 'Diamond Difference'
                for j = 1: length(Angle_positive) % use different pages to store the information for different angles
                    L_plus_T_positive(:, :, j) = Diagonal_ones;
                    L_plus_T_negative(:, :, j) = Diagonal_ones;
                    for i = 2: length(Spatial_Material_Parameter.Thickness)
                        L_plus_T_positive(i ,i-1, j) = - (1-tau_z(i, j)/2)/(1+tau_z(i, j)/2); 
                        L_plus_T_negative(i ,i-1, j) = - (1+tau_z(i, j)/2)/(1-tau_z(i, j)/2); 
                    end
                end
            otherwise
                disp('Spatial Discretization not available')
         end
        
        [Scalar_flux, First_angular_moment, existing_angular_flux_positive, existing_angular_flux_negative, psi_positive, psi_negative] = iterative_solver_1(psi_boundary, L_plus_T_positive, L_plus_T_negative, Spatial_Material_Parameter, Angle_positive, Weight_positive, Angle_negative, Weight_negative, tau_z, Spatial_discretization_option);
        
    case 2
        
        tau_reverse = flip(tau_z);
        switch Spatial_discretization_option
            case 'Step Characteristics'
                for j = 1: length(Angle_positive) % use different pages to store the information for different angles
                    L_plus_T_positive(:, :, j) = Diagonal_ones;
                    L_plus_T_negative(:, :, j) = Diagonal_ones;
                    for i = 2: length(Spatial_Material_Parameter.Thickness)
                        L_plus_T_positive(i ,i-1, j) = - exp(-tau_z(i, j));
                        L_plus_T_negative(i ,i-1, j) = - exp(-tau_reverse(i, j)); 
                    end
                end
            case 'Diamond Difference'
                for j = 1: length(Angle_positive) % use different pages to store the information for different angles
                    L_plus_T_positive(:, :, j) = Diagonal_ones;
                    L_plus_T_negative(:, :, j) = Diagonal_ones;
                    for i = 2: length(Spatial_Material_Parameter.Thickness)
                        L_plus_T_positive(i ,i-1, j) = - (1-tau_z(i, j)/2)/(1+tau_z(i, j)/2); % know psi_0, unknown, [psi_1; ... psi_N], omega>0
                        L_plus_T_negative(i ,i-1, j) = - (1-tau_reverse(i, j)/2)/(1+tau_reverse(i, j)/2); % know psi_N, unknown, [psi_N-1; ... psi_0], omega<0
                    end
                end
            otherwise
                disp('Spatial Discretization not available')
         end
        
        [Scalar_flux, First_angular_moment,  existing_angular_flux_positive, existing_angular_flux_negative, psi_positive, psi_negative] = iterative_solver_2(psi_boundary, L_plus_T_positive, L_plus_T_negative, Spatial_Material_Parameter, Angle_positive, Weight_positive, Angle_negative, Weight_negative, tau_z, Spatial_discretization_option);
        % existing_angular_flux_positive = psi_N(omega > 0)
        % existing_angular_flux_negative = psi_N(omega < 0)
    otherwise
        disp('unexpected boundary flag!')
end

end