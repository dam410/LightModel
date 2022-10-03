% Function to calculate source position using mirror
function [S,err,inters,normals] = calibrate_source(list_img,cameraParams,tagSize)
	% Right now use the default calibration
	intrinsics = cameraParams.Intrinsics;
	K = [intrinsics.FocalLength(1),0,intrinsics.PrincipalPoint(1);0,intrinsics.FocalLength(2),intrinsics.PrincipalPoint(2);0,0,1];
	invK = inv(K);
	rays = [];
	planes = [];
	nb_img = length(list_img);
	L_cell = cell(1,nb_img);
	for i=1:length(list_img)	
		I = imread(list_img{i});
		I = undistortImage(I,intrinsics,"OutputView","same");
		[planes_pose,list_points] = detect_apriltag_p(I,K,tagSize,1);
		h = figure('Name','Select the source position');
		imshow(I);
		hold on;
		P_proj_marqueurs = K*list_points{1};
		p_marqueurs = P_proj_marqueurs(1:2,:)./P_proj_marqueurs(3,:);
		Nd = planes_pose{1};
		N = Nd(1:3);
		d = Nd(4);
		plot(p_marqueurs(1,:),p_marqueurs(2,:),'-r');
		roi = drawpoint();
		pause;
		hold on;
		ray = invK*[transpose(roi.Position);1];
		lambda = -d/dot(N,ray);
		inter = lambda*ray;
		normal = inter-0.1*N;
		dir = inter-2*dot(N,inter)*N;
		dir = dir/norm(dir);
		pt_dir = inter+0.1*dir;
		proj_normal = K*normal;
		proj_normal = proj_normal(1:2)/proj_normal(3);
		proj_pt_dir = K*pt_dir;
		proj_pt_dir = proj_pt_dir(1:2)/proj_pt_dir(3);
		rays = [rays,invK*[transpose(roi.Position);1]];
		planes = [planes,[N;d]];
		plot([roi.Position(1),proj_normal(1)],[roi.Position(2),proj_normal(2)],'-b');
		plot([roi.Position(1),proj_pt_dir(1)],[roi.Position(2),proj_pt_dir(2)],'-r');
	end

	% First approximation Calculate the point of intersection between rays and plane
	lambdas = -planes(4,:)./dot(rays,planes(1:3,:));
	inters = lambdas.*rays;
	% Calculate the vector into reflected direction
	dirs = 2*inters-2*dot(inters,planes(1:3,:)).*planes(1:3,:);
	normals = -planes(1:3,:);
	normals_pts = inters+normals;
	figure('Name','Visualize all the rays with their plane intersection');
	plot3([zeros(1,length(list_img));inters(1,:)],...
		[zeros(1,length(list_img));inters(2,:)],...
		[zeros(1,length(list_img));inters(3,:)],...
		'-k');
	hold on;
	plot3([inters(1,:);dirs(1,:)],...
                [inters(2,:);dirs(2,:)],...
                [inters(3,:);dirs(3,:)],...
                '-r');
	plot3([inters(1,:);normals_pts(1,:)],...
                [inters(2,:);normals_pts(2,:)],...
                [inters(3,:);normals_pts(3,:)],...
                '-b');
	axis equal;
	[L_cell] = lines_2p_to_plucker(inters,dirs);
	[S,err] = closest_point_lines(L_cell)
end
