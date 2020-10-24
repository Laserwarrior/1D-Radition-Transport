function [Scalar_flux, First_angular_moment, existing_angular_flux_positive, existing_angular_flux_negative, psi_positive, psi_negative] = iterative_solver_2(psi_boundary, L_plus_T_positive, L_plus_T_negative, Spatial_Material_Parameter, Angle_positive, Weight_positive, Angle_negative, Weight_negative, tau_z, Spatial_discretization_option)

% both directions are defined on the left boundary

error_z_0 = 1; % L2 error initilized. This is not compared with the true solution or finest mesh solution, but with different iterations. 
error_z_1 = 1;

[phi_z_0, phi_z_1] = Initial_gaess_phi_z_i_2(psi_boundary, L_plus_T_positive, L_plus_T_negative, Spatial_Material_Parameter, Angle_positive, Weight_positive, Angle_negative, Weight_negative, tau_z, Spatial_discretization_option); % Based on no scattering. z represents the index for the material length = length(Spatial_Material_Parameter.Thickness), i represent ith order of angular moment for the zeroth spatial moment

% works before sigma_s < 1.7
% phi_z_0 = 0.5 + 0.5*ones(length(Spatial_Material_Parameter.Thickness), 1).*exp(3/5.*Spatial_Material_Parameter.Material_right_boundary);
% phi_z_1 = 0.5 + 0.5*ones(length(Spatial_Material_Parameter.Thickness), 1).*exp(3/5.*Spatial_Material_Parameter.Material_right_boundary);

% phi_z_0 = [2.9426
%     2.9999
%     3.0564
%     3.1120
%     3.1666
%     3.2203
%     3.2731
%     3.3249
%     3.3757
%     3.4255
%     3.4742
%     3.5220
%     3.5687
%     3.6144
%     3.6590
%     3.7025
%     3.7449
%     3.7862
%     3.8264
%     3.8654
%     3.9033
%     3.9401
%     3.9757
%     4.0101
%     4.0433
%     4.0753
%     4.1061
%     4.1357
%     4.1641
%     4.1913
%     4.2172
%     4.2418
%     4.2652
%     4.2874
%     4.3082
%     4.3278
%     4.3461
%     4.3632
%     4.3789
%     4.3933
%     4.4065
%     4.4183
%     4.4288
%     4.4381
%     4.4460
%     4.4526
%     4.4578
%     4.4618
%     4.4644
%     4.4657
%     4.4657
%     4.4644
%     4.4617
%     4.4577
%     4.4524
%     4.4458
%     4.4379
%     4.4287
%     4.4181
%     4.4062
%     4.3931
%     4.3786
%     4.3628
%     4.3458
%     4.3275
%     4.3078
%     4.2869
%     4.2648
%     4.2413
%     4.2167
%     4.1907
%     4.1636
%     4.1352
%     4.1055
%     4.0747
%     4.0426
%     4.0094
%     3.9750
%     3.9394
%     3.9026
%     3.8647
%     3.8256
%     3.7854
%     3.7441
%     3.7016
%     3.6581
%     3.6135
%     3.5678
%     3.5211
%     3.4733
%     3.4245
%     3.3747
%     3.3238
%     3.2720
%     3.2193
%     3.1655
%     3.1109
%     3.0553
%     2.9988
%     2.9414]*10^3;
% 
% phi_z_1 = [-2.8897
%    -2.8458
%    -2.8010
%    -2.7553
%    -2.7088
%    -2.6616
%    -2.6135
%    -2.5647
%    -2.5151
%    -2.4648
%    -2.4137
%    -2.3620
%    -2.3095
%    -2.2564
%    -2.2025
%    -2.1481
%    -2.0930
%    -2.0372
%    -1.9809
%    -1.9240
%    -1.8665
%    -1.8085
%    -1.7499
%    -1.6908
%    -1.6312
%    -1.5711
%    -1.5106
%    -1.4496
%    -1.3882
%    -1.3263
%    -1.2641
%    -1.2015
%    -1.1386
%    -1.0753
%    -1.0117
%    -0.9478
%    -0.8836
%    -0.8191
%    -0.7545
%    -0.6895
%    -0.6244
%    -0.5591
%    -0.4937
%    -0.4280
%    -0.3623
%    -0.2965
%    -0.2305
%    -0.1645
%    -0.0985
%    -0.0324
%     0.0337
%     0.0998
%     0.1658
%     0.2318
%     0.2978
%     0.3636
%     0.4293
%     0.4950
%     0.5604
%     0.6257
%     0.6908
%     0.7557
%     0.8204
%     0.8849
%     0.9490
%     1.0129
%     1.0765
%     1.1398
%     1.2028
%     1.2654
%     1.3276
%     1.3894
%     1.4508
%     1.5118
%     1.5723
%     1.6324
%     1.6920
%     1.7510
%     1.8096
%     1.8676
%     1.9251
%     1.9820
%     2.0383
%     2.0940
%     2.1491
%     2.2036
%     2.2574
%     2.3106
%     2.3630
%     2.4148
%     2.4658
%     2.5161
%     2.5657
%     2.6145
%     2.6625
%     2.7098
%     2.7562
%     2.8018
%     2.8466
%     2.8906]*10^3;

phi_z_0_positive = zeros(length(Spatial_Material_Parameter.Thickness), length(Angle_positive));
phi_z_1_positive = phi_z_0_positive;
phi_z_0_negative = phi_z_0_positive;
phi_z_1_negative = phi_z_0_positive;
psi_positive = phi_z_0_positive;
psi_negative = phi_z_0_positive;

existing_angular_flux_positive = zeros(1, length(Angle_positive));
existing_angular_flux_negative = existing_angular_flux_positive;

tau_reverse = flip(tau_z);

while error_z_0 > 10^(-6) && error_z_1 > 10^(-6)
    
    for i = 1:length(Angle_positive)
        %% positive directions
        [phi_z_0_positive(:, i), phi_z_1_positive(:, i), existing_angular_flux_positive(i), psi_positive(:, i)] = sovle_matrix_known_phi(phi_z_0,       phi_z_1,       Angle_positive(i), Spatial_Material_Parameter,       tau_z(:, i),       psi_boundary(i,1), L_plus_T_positive(:, :, i), 1, Weight_positive(i), Spatial_discretization_option);
        % for the existing_angular_flux_negative, it is actually not existing
        % psi_positive = [psi_1, psi_2,..., psi_N]
        % but incoming. So, it might not be considered when presenting the data
        %% negative directions
        [phi_z_0_negative(:, i), phi_z_1_negative(:, i), existing_angular_flux_negative(i), psi_negative(:, i)] = sovle_matrix_known_phi(flip(phi_z_0), flip(phi_z_1), Angle_negative(i), flip(Spatial_Material_Parameter), tau_reverse(:, i), psi_boundary(i,2), L_plus_T_negative(:, :, i), 1, Weight_negative(i), Spatial_discretization_option);
        % psi_negative = [psi_N-1, psi_N-2,..., psi_0]
    end
    
    phi_z_0_new = sum(phi_z_0_positive, 2) + sum(flip(phi_z_0_negative), 2);
    phi_z_1_new = sum(phi_z_1_positive, 2) + sum(flip(phi_z_1_negative), 2);
    
    error_z_0 = sqrt( sum( (phi_z_0_new - phi_z_0).^2 ) );
    error_z_1 = sqrt( sum( (phi_z_1_new - phi_z_1).^2 ) );
    
    phi_z_0 = phi_z_0_new;
    phi_z_1 = phi_z_1_new;
    
end

Scalar_flux = phi_z_0; % zeroth spatial moment and zeroth angular moment
First_angular_moment = phi_z_1; % zeroth spatial moment and first angular moment


end