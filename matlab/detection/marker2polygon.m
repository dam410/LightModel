function [mesh_] = marker2polygon(cell_R_i,cell_t_i,tagSize,isboard,options)
%(cell_R_i,cell_t_i,tagSize,isboard)
	arguments
		cell_R_i (1,:) cell
		cell_t_i (1,:) cell
		tagSize (1,1) double
		isboard (1,1) double
		options.specific_coordinates (2,:) double
	end
	% For each marker we calculate 3 points
	nb_marqueur = length(cell_R_i);
	P = [];
	mesh_ = struct();
	mesh_.pose = eye(4);
	mesh_.type = 'mesh';
	data = struct();
	P_marqueur = [];
	N_marqueurs = [];
	for j=1:nb_marqueur
		R = transpose(cell_R_i{j});
		t = cell_t_i{j};
		N_marqueurs = [N_marqueurs,-R(:,3)];
		P = [P,2*tagSize*R(:,1)+t,2*tagSize*R(:,2)+t,t];
		%P = [P,t];
		P_marqueur = [P_marqueur,t];
	end
	% Estimate the plane with points + help from orthogonal vectors
	%	using marker coordinate systems
	%N_marqueurs
	%N = median(N_marqueurs,2);	
	%N = N/norm(N)
	%d = -mean(transpose(N)*P_marqueur);
	%Nd = [N;d];
	Nd = estimate_plane(P);
	% Setting default material
	material = struct();
	material.ambient = 0;
	material.diffuse_color = [1;1;1];
	material.diffuse_intensity = 0.7;
	material.diffuse_shader = 'LAMBERT';
	material.specular_color = [1;1;1];
	material.specular_intensity = 0;
	material.specular_shader = 'PHONG';
	data.materials = material;
	% Setting normals
	data.normals = ones(nb_marqueur,1)*transpose(Nd(1:3));
	% Projecting points in 2D and calculating the convex hull to extract a polygon
	R_plane = [null(Nd(1:3)*transpose(Nd(1:3))),Nd(1:3)];
	if det(R_plane)<0
		R_plane = R_plane(:,[2,1,3]);
	end
	p_marqueur = transpose(R_plane)*P_marqueur;
	p_marqueur = p_marqueur(1:2,:);
	k = convhull(transpose(p_marqueur));
	% If we use board we translate the point to get only the central mesh
	if isboard>0
		P_mid = mean(P_marqueur,2);
		dir_P_mid = (P_mid-P_marqueur)./sqrt(sum((P_mid-P_marqueur).^2));
		P_translated = P_marqueur+isboard*tagSize*dir_P_mid;
		data.vertices = transpose(P_translated);
		%figure('Name','Testing  the translation');
		%plot3(P_marqueur(1,:),P_marqueur(2,:),P_marqueur(3,:),'+r');
		%hold on;
		%plot3(P_translated(1,:),P_translated(2,:),P_translated(3,:),'+g');
	else
		data.vertices = transpose(P_marqueur);
	end
	polygon = struct();
	polygon.material_index = 0;
	polygon.normal = Nd(1:3);
	% From Matlab convention for starting index to C/Python
	polygon.vertices_index = k(1:(end-1))-1;

	% If we have the option to select a specific polygonal design
	% Remake the entire polygon
	if isfield(options,'specific_coordinates')
		% Specific coordinates are barycentric coordinates based on the markers corners
		% See file calculate_polgon_markers in general_function
		p_poly_coord = options.specific_coordinates;
		data.normals = transpose(p_poly_coord)*data.normals(1:2,:);
		data.vertices = transpose(p_poly_coord)*(data.vertices(1:2,:)-data.vertices(3,:))+...
			data.vertices(3,:);
		polygon.vertices_index = 0:1:(length(p_poly_coord)-1);

	end
	
	data.polygons = [polygon];
	mesh_.data = data;
end
