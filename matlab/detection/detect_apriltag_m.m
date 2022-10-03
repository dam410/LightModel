% Detect all the apriltag in one image and return it like a mesh
function [meshes_vector] = detect_apriltag_m(I,K,tagSize,nbPlan,isboard)

	if nargin<5
		isboard = 0;
	end

	filepath = 'temp.pgm';
	imwrite(I,filepath);
	[cell_R,cell_t,cell_id,cell_err] = detect_apriltag(filepath,K,tagSize);
	disp('Nb marker detected, detect_apriltag_m: line 11');
	length(cell_R)
	cell_R = cell_R(cell_err<1e-6);
	cell_t = cell_t(cell_err<1e-6);
	cell_id = cell_id(cell_err<1e-6);
	cell_err = cell_err(cell_err<1e-6);


	% Segmenter les marqueurs par plan puis réestimer à partir de tous les marqueurs et les compter
	% 	comme des polygones
	P = [];
	for i=1:length(cell_R)	
		R = transpose(cell_R{i});
		t = cell_t{i};
		N = R(:,3);
		d = -dot(N,t);
		P = [P,[N;d]];
	end
	% Creating output meshes_vector
	meshes_vector = [];
	if nbPlan > 1
		[idx] = kmeans(transpose(P),nbPlan);
		for n_plane = 1:nbPlan
			cell_R_i = cell_R(idx==n_plane);
			cell_t_i = cell_t(idx==n_plane);
			meshes_vector = [meshes_vector,marker2polygon(cell_R_i,cell_t_i,tagSize,isboard)];
		end
	else
			meshes_vector = [meshes_vector,marker2polygon(cell_R,cell_t,tagSize,isboard)];
	end
end
