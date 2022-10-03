function [isocontour_pts,levels,curve_params] = detect_isocontour(I,options)
	arguments
		I (:,:) double
		options.Polygon (:,:) double
		options.StandardDeviation (1,1) double = 3
		options.FilterSize (1,1) uint8  = 11
		options.Shrink (1,1) double = 3
		options.NbIsocontour (1,1) double = 2
		options.Fitting (1,1) string
		options.NbIsoArtificial (1,1) = 20
	end
	% Find the contour by using ccontour of Matlab with gaussian blur
	I_gaussian = imgaussfilt(I,options.StandardDeviation,'FilterSize',options.FilterSize);
	% Restrict the area to the polygon if provided
	if isfield(options,'Polygon')
		poly_2D = [options.Polygon,options.Polygon(:,1)];
		%	Calculating the mask
		[n,m] = size(I)
		BW = poly2mask(poly_2D(1,:),poly_2D(2,:),n,m);
		% 	Shrinking allow to remove inaccuracy on the blender projection
		%BW = bwmorph(BW,'shrink',options.Shrink);
		I_gaussian(~BW) = 0;
		% Another way to do it is to calculate a smaller polygon, do the normal
		% detection or masked then we remove points that are outside the smaller
		% polygon.
		%	Calculating the smaller polygon
		pt_barycenter = mean(poly_2D,2);
		direction_shrink = pt_barycenter-poly_2D;
		norm_direction_shrink = sqrt(sum(direction_shrink.^2));
		direction_shrink = direction_shrink./norm_direction_shrink;
		poly_2D_smaller = poly_2D+15*direction_shrink;
	end
	% Parameters to change the real number of isocontours considered
	NB_ISO_ADD = options.NbIsoArtificial;
	C = contourc(I_gaussian,options.NbIsocontour+NB_ISO_ADD);
	%C = contourc(I_gaussian);
	% Read the points from the contour found
	[points,levels] = contour_to_points(C);
	NB_ISO_ADD = length(points)-options.NbIsocontour;
	nb_pts = [];
	for i=1:length(points)
		nb_pts = [nb_pts,length(points{i})];
	end
	
	curve_params = {};
	isocontour_pts = {};
	for i = NB_ISO_ADD:(length(levels)-1)
		data_points = [transpose(points{i}(2,:)),transpose(points{i}(1,:))];
		if isfield(options,'Polygon')
			trans_data_points = transpose(data_points);
			in_poly = inpolygon(trans_data_points(2,:),trans_data_points(1,:),...
				poly_2D_smaller(1,:),poly_2D_smaller(2,:));
			data_points = transpose(trans_data_points(:,in_poly));
			%figure('Name','test');
			%plot(poly_2D_smaller(1,:),poly_2D_smaller(2,:),'-b');
			%hold on;
			%plot(trans_data_points(2,:),trans_data_points(1,:),'+r');
		end
		if isfield(options,'Fitting')
			switch options.Fitting
			case 'ellipse'
				[ellipse_param,~] = ellipseFromPoints(data_points);
				curve_params{i-NB_ISO_ADD+1} = ellipse_param;
			otherwise
				curve_params{i-NB_ISO_ADD+1} = [];
			end
		end
		isocontour_pts{i-NB_ISO_ADD+1} = data_points;
	end
	levels = levels(NB_ISO_ADD:(length(levels)-1));
end

