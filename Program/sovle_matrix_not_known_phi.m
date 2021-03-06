function [phi_z_0, phi_z_1] = sovle_matrix_not_known_phi(Angle_single_direction, Spatial_Material_Parameter, tau_z, psi_boundary, L_plus_T_single_direction, direction_flag, weight, Spatial_discretization_option)
     
    % used to initialize the phi_z_0 and phi_z_1
    % here, single direction means either left or right, mean only one
    % angle
 
    
    Q_single_direction = Spatial_Material_Parameter.Q_0 * ones(1,length(Angle_single_direction)) + Spatial_Material_Parameter.Q_1 * Angle_single_direction' ; % each collum represents different directions, each row represents spatial ordinate

    Q_star_single_direction = Q_single_direction; % each collum represents different directions, each row represents different positions
    
    switch Spatial_discretization_option
    case 'Step Characteristics'
        b_Q_star_single_direction = Q_star_single_direction./ Spatial_Material_Parameter.Sigma_t .* (1 - exp(-tau_z)); % still each collum represents different directions, each row represents different position
        b_Q_star_single_direction(1,:) = b_Q_star_single_direction(1,:) + psi_boundary' * exp(-tau_z(1)); % here, psi_boundary is just a column array
        
    case 'Diamond Difference'
        b_Q_star_single_direction = Q_star_single_direction./ Spatial_Material_Parameter.Sigma_t .* (1 - (1-tau_z/2)./(1+tau_z/2)); % still each collum represents different directions, each row represents different position
        b_Q_star_single_direction(1,:) = b_Q_star_single_direction(1,:) + psi_boundary' * (1-tau_z(1)/2)/(1+tau_z(1)/2); % here, psi_boundary is just a column array
    otherwise
        disp('Spatial Discretization not available')
    end

    psi_single_direction = L_plus_T_single_direction\b_Q_star_single_direction;
    psi_single_direction_index_minus_one = zeros(size(psi_single_direction));
    psi_single_direction_index_minus_one(2:end,1:end) = psi_single_direction(1:end-1, 1:end);
    psi_single_direction_index_minus_one(1,1:end) = psi_boundary';
    
    phi_z_0 = sum( ((psi_single_direction_index_minus_one - psi_single_direction)./(direction_flag*tau_z) + Q_star_single_direction./ Spatial_Material_Parameter.Sigma_t).*weight', 2);
    phi_z_1 = sum( ((psi_single_direction_index_minus_one - psi_single_direction)./(direction_flag*tau_z) + Q_star_single_direction./ Spatial_Material_Parameter.Sigma_t).*Angle_single_direction'.*weight', 2);


end