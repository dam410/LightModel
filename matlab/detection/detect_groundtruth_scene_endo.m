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
		options.Polygon (1,1) double = 0
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
		poly_coord = [-0.1204,-0.0983, 0.1644, 0.1765, 0.8067, 0.7947, 1.0787, 1.0592, 0.8241, 0.8071, 0.1646, 0.1588, -0.1204;...
		0.9743, 0.2521,-0.0086,-0.2697,-0.8910,-0.6378,-0.8947,-0.2057, 0.0137, 0.2965, 0.9277, 0.7074, 0.9743];
		if options.Polygon
			meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,0,'specific_coordinates',poly_coord)];
		else
			meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
		end
	end


	% Sticker black 2 -> cell_id = [26 28 29 4]
	[Lia,Locb] = ismember([26 28 29 4],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		poly_coord = [0.7904,0.1501,0.1652,-0.1252,-0.1237,0.1310,0.1573,0.7807,0.7773,1.0630,1.0593,0.7911, 0.7904;...
		0.2937, 0.9368, 0.6897, 0.9955, 0.2645, 0.0167,-0.2512,-0.8715,-0.6607,-0.9234,-0.2099, 0.0460,0.2937];
		if options.Polygon
			meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,0,'specific_coordinates',poly_coord)];
		else
			meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
		end
	end


	% Sticker brown -> cell_id = [20 21 10 25]
	[Lia,Locb] = ismember([20 21 10 25],cell_id);
	Locb = Locb(find(Locb>0));
	if length(find(Locb>0))>0
		poly_coord = [1.0874, 1.0831, 0.7803, 0.7778, 0.1890, 0.2250,-0.1194,-0.1086, 0.1850, 0.2151, 0.7960, 0.7911, 1.0874;...
		-0.9033,-0.2366, 0.0642, 0.3025, 0.8672, 0.6613, 0.9699, 0.2346,-0.0319,-0.2941,-0.8425,-0.6316, -0.9033];
		if options.Polygon
			meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,0,'specific_coordinates',poly_coord)];
		else
			meshes_vector = [meshes_vector,marker2polygon(cell_R(Locb),cell_t(Locb),tagSize,1.7)];
		end
	end

end

