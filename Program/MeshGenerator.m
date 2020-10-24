function Spatial_Material_Parameter = MeshGenerator(Material_distribution, Mesh_error_index, i, Macroscopic_Cross_Sections)

Material_Mesh = zeros(length(Material_distribution.Material),1);

Material_error_analysis_mesh = Material_distribution.Mesh{Mesh_error_index};
Material_Mesh(Mesh_error_index) = Material_error_analysis_mesh(i);

if length(Material_distribution.Material) > 1
    
    switch Mesh_error_index
    
        case 1
            for i = 2:length(Material_distribution.Material)
                Material_Mesh(i) = Material_distribution.Mesh{i};
            end
        
        case length(Material_distribution.Material)
            for i = 1:length(Material_distribution.Material)-1
                Material_Mesh(i) = Material_distribution.Mesh{i};
            end
        
        otherwise
            for i = 1:Mesh_error_index-1
                Material_Mesh(i) = Material_distribution.Mesh{i};
            end
            
            for i = Mesh_error_index+1:length(Material_distribution.Material)
                Material_Mesh(i) = Material_distribution.Mesh{i};
            end
    end
    
end

Material = cell( sum(Material_Mesh) ,1);
Thickness = zeros( sum(Material_Mesh) ,1 );
Material_Mesh_right_boundary = cumsum(Material_Mesh);

Material(1: Material_Mesh(1)) = Material_distribution.Material(1);
Thickness(1: Material_Mesh(1)) = Material_distribution.Thickness(1) / Material_Mesh(1);

if length(Material_Mesh) > 1
    
    for i = 2:length(Material_Mesh)
        
        Material( Material_Mesh_right_boundary(i-1)+1: Material_Mesh_right_boundary(i)) = Material_distribution.Material(i);
        Thickness( Material_Mesh_right_boundary(i-1)+1: Material_Mesh_right_boundary(i)) = ones(Material_Mesh(i),1) * Material_distribution.Thickness(i) / Material_Mesh(i);
    
    end
    
end

Spatial_Material_Parameter = Macroscopic_Cross_Sections(Material,:);
Spatial_Material_Parameter.Thickness = Thickness;
Spatial_Material_Parameter.Material_right_boundary = cumsum(Spatial_Material_Parameter.Thickness);

end