function [Angle_discretization_quadrature_set_name, Angle_discretization_quadrature_set_order] = Anlge_discretization_method(Anlge_discretization_option)

colon_position = strfind(Anlge_discretization_option,':');
Angle_discretization_quadrature_set_name = Anlge_discretization_option(1:colon_position - 1);
Angle_discretization_quadrature_set_order = Anlge_discretization_option(colon_position+2:end);


end