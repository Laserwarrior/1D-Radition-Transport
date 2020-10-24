function [psi_boundary, boundary_flag] = psi_boundary_fun(Surface_flag, Surface_Sources, Angle_positive, Angle_negative, Angle)

if Surface_flag == 1 % left boundary is not vacuum. In the setting, the right will always be vacuum
    psi_0_fun = str2func( char ( Surface_Sources.(2) ) );
    Range_of_Omega_z = char( Surface_Sources.(1) );
    
    if strcmp(Range_of_Omega_z, 'None') % Left boundary covers omega_z from -1 to 1. Don't need to consider the right boundary.
        boundary_flag = 1;
        if strcmp(Surface_Sources.Row{1}, 'Discrete') 
            psi_0 = arrayfun(psi_0_fun, sym(Angle)); % discrete (because of dirac function, still only need to do the sweep from left to right)
            psi_0_positive = psi_0(2:2:end);
            psi_0_negative = psi_0(1:2:end-1);
            psi_boundary = [psi_0_positive, psi_0_negative];
            psi_boundary = double(psi_boundary);          
        else                                          
            psi_0 = arrayfun(psi_0_fun, Angle);                  % linear (sweep two directions based on the left boundary, since both directions are defined on the left boundary)
            psi_0_positive = psi_0(2:2:end);
            psi_0_negative = psi_0(1:2:end-1);
            psi_boundary = [psi_0_positive, psi_0_negative];
        end
        
    else  % Isotropic. Left boundary doesn't covers omega_z from -1 to 1. Need to consider the right boundary.
        psi_0_positive = arrayfun(psi_0_fun, Angle_positive); % Isotropic; Cover omega_z from 0 to 1
        psi_N_negative = zeros(length(Angle_positive), 1); % In the problem, right is vacuum. So the BC on the right going to left should be 0 Cover omega_z from -1 to 0
        psi_boundary = [psi_0_positive, psi_N_negative];
        boundary_flag = 2;
    end
    
else
    psi_0_positive = zeros( length(Angle)/2, 1 ); % on the left going to right is 0. Omega_z from  0 to 1
    psi_N_negative = zeros( length(Angle)/2, 1 ); % on the right going to left is 0. Omega_z from -1 to 1
    psi_boundary = [psi_0_positive, psi_N_negative];
    boundary_flag = 2;
end


end