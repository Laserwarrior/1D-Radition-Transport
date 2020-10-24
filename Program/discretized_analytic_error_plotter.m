function discretized_analytic_error_plotter(Scalar_flux, phi_analytic, psi_positive, Analytic_psi, Spatial_Material_Parameter)

x_coordinate = Spatial_Material_Parameter.Material_right_boundary;
figure;
plot(x_coordinate, Analytic_psi,'LineStyle','-','Color','g','LineWidth',6); hold on

error_psi = abs(Analytic_psi - psi_positive);
errorbar(x_coordinate, psi_positive, error_psi, 'or'); 

figure;
plot(x_coordinate, phi_analytic, 'LineStyle','-','Color','g','LineWidth',6); hold on

error_phi = abs(phi_analytic - Scalar_flux);
errorbar(x_coordinate, Scalar_flux, error_phi, 'or'); 



end