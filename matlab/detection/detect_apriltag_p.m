% Detect all the apriltag in one image and return only the pose of the plane
% 	and the position of the tag detected
function [planes_pose,list_points] = detect_apriltag_p(I,K,tagSize,nbPlan)

	filepath = 'temp.pgm';
	imwrite(I,filepath);
	[cell_R,cell_t,cell_id,cell_err] = detect_apriltag(filepath,K,tagSize);
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
		planes_pose = cell(1,nbPlan);
		for n_plane = 1:nbPlan
			% Reject outliers by checking more than
			N_med = median(P(1:3,idx==n_plane),2);
			N_med = N_med/norm(N_med);
			check_angle = acos(dot(repmat(N_med,1,size(P,2)),P(1:3,:)))*180/pi < 10
			cell_R_i = cell_R(idx==n_plane & check_angle);
			cell_t_i = cell_t(idx==n_plane & check_angle);
			[Nd,P_marqueur] = marker2plane_pose(cell_R_i,cell_t_i,tagSize);
			
		end
	else
			% Reject outliers by checking more than
			N_med = median(P(1:3,:),2);
			N_med = N_med/norm(N_med);
			check_angle = acos(dot(repmat(N_med,1,size(P,2)),P(1:3,:)))*180/pi < 10
			[Nd,P_marqueur] = marker2plane_pose(cell_R(check_angle),cell_t(check_angle),tagSize);
			planes_pose = {Nd};
			list_points = {P_marqueur};
	end
end

function [Nd,P_marqueur] = marker2plane_pose(cell_R_i,cell_t_i,tagSize)
	% For each marker we calculate 3 points
	nb_marqueur = length(cell_R_i);
	P = [];
	P_marqueur = [];
	N_marqueurs = [];
	for j=1:nb_marqueur
		R = transpose(cell_R_i{j});
		t = cell_t_i{j};
		N_marqueurs = [N_marqueurs,-R(:,3)];
		P = [P,2*tagSize*R(:,1)+t,2*tagSize*R(:,2)+t,-2*tagSize*R(:,1)+t,-2*tagSize*R(:,2)+t,t];
		%P = [P,t];
		P_marqueur = [P_marqueur,t];
	end
	% Estimate the plane with points + help from orthogonal vectors
	%	using marker coordinate systems
	Nd = estimate_plane(P);
end
