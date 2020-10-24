function [Scalar_flux, First_angular_moment, existing_angular_flux_positive, Angle_positive, psi_positive, Angle_negative, psi_negative, psi_boundary] = Solver(Surface_Sources, Surface_flag, Spatial_Material_Parameter, Anlge_discretization_option, Spatial_discretization_option, Solving_method)

switch Solving_method
    case 'Transport_Sweep'
        [Scalar_flux, First_angular_moment, existing_angular_flux_positive, Angle_positive, psi_positive, Angle_negative, psi_negative, psi_boundary] = Transport_Sweep_solver(Surface_Sources, Surface_flag, Spatial_Material_Parameter, Anlge_discretization_option, Spatial_discretization_option) ;
        
    otherwise
        disp('Solving Method Not Available')
end

end