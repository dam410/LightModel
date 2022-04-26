% Calculate the pose of each scene element using markers
%	Right now we can support this elements:
%
% Sticker Large -> cell_id = [2 22 26 28]
% Canson Small -> cell_id = |8 13 18 23]
% Canson Large -> cell_id = |7 12 17 27]
% Candle board -> cell_id = |0 1 5 6]

function [meshes_vector,S] = detect_groundtruth_scene(I_undistorted,K,options)
	arguments
                I_undistorted (:,:,:) uint8
                K (3,3) double
                options.CandleHeight (1,1) double = 0.095
                options.TagSize (1,1) double = 0.03
        end
	tagSize = options.TagSize;
	filepath = 'temp.pgm';
        imwrite(I_undistorted,filepath);
	K = K*diag([1,1,1]);
        [cell_R,cell_t,cell_id,cell_err] = detect_apriltag(filepath,K,tagSize);
	% Find the scene elements using their id
	meshes_vector = [];	
	S = [0;0;0];
	% Sticker Large -> cell_id = [2 22 26 28]
	[Lia,Locb] = ismember([2 22 26 28],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
	end

	% Canson Small -> cell_id = [8 13 18 23]
	[Lia,Locb] = ismember([8 13 18 23],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
	end

	% Canson Large -> cell_id = [7 12 17 27]
	[Lia,Locb] = ismember([7 12 17 27],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
	end

	% Candle board -> cell_id = [0 1 5 6]
	[Lia,Locb] = ismember([0 1 5 6],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		S = marker2source(cell_R(Locb),cell_t(Locb),tagSize,options.CandleHeight);
	end

	%[h,l,~] = size(I_undistorted);
	%[I_test,~,~,~,P_camera] = render_shading_isocontour(l,h,...
        %                'Surface','Meshes',...
        %                'Meshes',meshes_vector,...
        %                'Scattering','Phong',...
        %                'LightParameters',[S;0.6],...
        %                'CameraIntrinsic',K,...
        %                'UseImageTransformation',0);
	%I_test = permute(I_test,[2 1 3]);

	%figure('Name','Simulation of the image according model');
	%imshow(I_test);

end

function [S] = marker2source(cell_R,cell_t,tagSize,height)
	% Similarly to the board we calculate the plane and its normal
	P = [];
	P_markers = [];
	for j=1:length(cell_R)
                R = transpose(cell_R{j});
                t = cell_t{j};
		P_markers = [P_markers,t];
                P = [P,2*tagSize*R(:,1)+t,2*tagSize*R(:,2)+t,t];
                %P = [P,t];
        end
	Nd = estimate_plane(P);
	% If all the markers are found, we can easily take the barycenter for the 2D coordinate of the
	%	source and add
	if length(cell_R) == 4
		S = height*Nd(1:3)+mean(P_markers,2);
	elseif length(cell_R) == 3
		P_3 = P_markers;
		% Find the right-angled in 3 points
		vect_3 = P_3(:,1:3) - P_3(:,[2,3,1]);
		vect_3 = vect_3./sqrt(sum(vect_3.^2));
		dot_vect_3 = acos(abs(dot(vect_3(:,1:3),vect_3(:,[3,1,2]))));
		[~,index_45_45_90] = sort(dot_vect_3,'ascend');
		% Find the middle point of the two other points
		S = height*Nd(1:3) + 0.5*P_3(:,index_45_45_90(1)) + 0.5*P_3(:,index_45_45_90(2));
	end
end

