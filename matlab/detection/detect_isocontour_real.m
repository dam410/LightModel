function [isocontour_pts,levels,curve_params] = detect_isocontour_real(I_rgb,polys_2D,nbIso,display)


	if size(I_rgb,3)==3
		I_double = double(rgb2gray(uint8(I_rgb)));
	else
		I_double = double(I_rgb);
	end
	I_double = imgaussfilt(I_double,3,'FilterSize',11);
	% Filter to remove strong noise
	I_double = wiener2(I_double,[51,51]);
	I_double = imgaussfilt(I_double,3,'FilterSize',11);
	[n,m] = size(I_double);


	% Intensity in percentage that will be used for isocontour and ellipse detection
	%p_tested = [10:17]/20;
	interval = 0:1/nbIso:1;
	p_tested = 14/20*(interval(1:(end-1))-mean(interval(1:(end-1))))+13/20;

	curve_params = cell(length(polys_2D),length(p_tested));
	isocontour_pts =  cell(length(polys_2D),length(p_tested));
	levels = cell(length(polys_2D),length(p_tested));
	
	% Extract orthogonal lines through isocontour
	length(polys_2D)
	for i_poly = 1:length(polys_2D)
		I_current = I_double;
	
		% Remove the image part outside the polygon
		poly_2D = [polys_2D{i_poly},polys_2D{i_poly}(:,1)];
		pt_barycenter = mean(poly_2D,2);
                direction_shrink = pt_barycenter-poly_2D;
                norm_direction_shrink = sqrt(sum(direction_shrink.^2));
                direction_shrink = direction_shrink./norm_direction_shrink;
                poly_2D_smaller = poly_2D+15*direction_shrink;
		BW = poly2mask(poly_2D_smaller(1,:),poly_2D_smaller(2,:),n,m);
		I_current(~BW) = 0;
		
		% Calculate the min and max intensity to adjust the isocontour detection
		min_I = min(I_current(I_current(:)>1));
		max_I = max(I_current(:));
		
		for i=1:length(p_tested)
			p_1 = p_tested(i)-1/600;
			p_2 = p_tested(i)+1/600;
			I_p_1 = p_1*(max_I-min_I)+min_I;
			I_p_2 = p_2*(max_I-min_I)+min_I;

			index_full = find((I_current > I_p_1 & I_current < I_p_2));
			[ind_x,ind_y] = ind2sub(size(I_current),index_full);
			data_points = [ind_x,ind_y];
			[ellipse_param,cov_ell_param] = ellipseFromPoints(data_points);
			isocontour_pts{i_poly,i} = data_points;
			curve_params{i_poly,i} = ellipse_param;
			levels{i_poly,i} = (I_p_1+I_p_2)/2;
		end
	end
	% Adding some display in case of failing
	if nargin > 3 && strcmp(display,'on')
		figure('Name','Show the image and detected isocontours');
		imshow(I_double/255);
		hold on;
		for i_poly = 1:length(polys_2D)
			for i = 1:nbIso
				plot(isocontour_pts{i_poly,i}(:,2),isocontour_pts{i_poly,i}(:,1),'+r');
			end
		end
	end

end

