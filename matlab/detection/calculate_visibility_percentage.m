function [eccentricity,visibility] = calculate_visibilty_percentage(ell,pts)
	% Calculate the percentage in angles otherwise it will be complicated for nothing	
	pts_centered = pts - [ell(1),ell(2)];
	T = [cos(ell(5)),sin(ell(5));-sin(ell(5)),cos(ell(5))]*transpose(pts_centered);
	T(1,:) = T(1,:)./ell(3);
	T(2,:) = T(2,:)./ell(4);
	ang_pts = atan2(T(2,:),T(1,:))+2*pi;
	ang_pts = sort([ang_pts-min(ang_pts),2*pi]);
	if length(ang_pts)>1
		visibility = (2*pi-max(ang_pts(2:end)-ang_pts(1:(end-1))))/(2*pi)
	else
		visibility = 0;
	end
	eccentricity = sqrt(1-ell(4)^2/ell(3)^2);
end
