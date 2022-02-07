function [Xcs,Ns] = estimate_plane_pose(data)
	n_p = length(data.homography);
	Xcs = cell(1,n_p);
	Ns = cell(1,n_p);
	for i_p = 1:n_p
		% Calculate the normal and the center using the homography
		nb_ambig = length(data.homography{i_p});
		N = zeros(3,nb_ambig);
		Xc = zeros(2,nb_ambig);
		for i_ambig = 1:nb_ambig
			H = data.homography{i_p}{i_ambig};
			invH = inv(H);
			N(:,i_ambig) = transpose(data.K)*transpose(H)*[0;0;1];
			N(:,i_ambig) = N(:,i_ambig)/norm(N(:,i_ambig));
			Xc_3D = invH*[0;0;1];
			Xc(:,i_ambig) = Xc_3D([1,2])/Xc_3D(3);
		end
		Xcs{i_p} = Xc;
		Ns{i_p} = N;
	end
end
