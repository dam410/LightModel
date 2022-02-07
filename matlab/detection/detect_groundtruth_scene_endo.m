% Calculate the pose of each scene element using markers
%	Right now we can support this elements:
%
% Sticker Large -> cell_id = [2 22 26 28]
% Canson Small -> cell_id = |8 13 18 23]
% Canson Large -> cell_id = |7 12 17 27]
% Candle board -> cell_id = |0 1 5 6]

function [meshes_vector,S] = detect_groundtruth_scene_endo(I_undistorted,K,options)
	arguments
                I_undistorted (:,:,:) uint8
                K (3,3) double
                options.TagSize (1,1) double = 0.008
        end
	tagSize = options.TagSize;
	filepath = 'temp.pgm';
        imwrite(I_undistorted,filepath);
	K = K*diag([1,1,1]);
        [cell_R,cell_t,cell_id,cell_err] = detect_apriltag(filepath,K,tagSize);
	% Find the scene elements using their id
	meshes_vector = [];	
	S = [0;0;0];
	% Sticker black 1 -> cell_id = [0 5 27 16]
	[Lia,Locb] = ismember([0 5 27 16],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
	end

	% Sticker black 2 -> cell_id = [26 28 29 4]
	[Lia,Locb] = ismember([26 28 29 4],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
	end

	% Sticker brown -> cell_id = [20 21 10 25]
	[Lia,Locb] = ismember([20 21 10 25],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
	end
end

