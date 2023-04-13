function leaflet = import_valve_from_csv(leaflet, file_name, indices_global, directory)

if exist('directory', 'var')
    file_with_path = fullfile(directory, file_name); 
else 
    file_with_path = file_name; 
end 

vertices = csvread(file_with_path); 

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 

for j = 1:j_max
    for k = 1:k_max
        vertex_idx = indices_global(j,k) + 1; 
        leaflet.X(:,j,k) = vertices(vertex_idx,:)' ; 
    end 
end 




