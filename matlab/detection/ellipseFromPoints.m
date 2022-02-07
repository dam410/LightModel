% This function compute an ellipse from points using algorithm proposed by Szpak.

function [ellipse_param,cov_ell_param,ellipse_norm,cov_ell_norm,T_param2norm] = ellipseFromPoints(data_points,pts_var)

	nb_points = size(data_points,1);
	covList = mat2cell(repmat(eye(2),1,nb_points),2,2.*(ones(1,nb_points)));
	if  nb_points == 0
		ellipse_param = zeros(1,5);
		cov_ell_param = zeros(5,5);
		return	
	end
	[theta_fastguaranteed] = fast_guaranteed_ellipse_estimate(data_points);
	if nargin < 2
		pts_var = estimateNoiseLevel(theta_fastguaranteed,...
		  data_points, covList);
	end
	for iPt = 1:nb_points
		covList{iPt} = pts_var*covList{iPt};
	end


	% If required we can also compute the ellipse in normalized 6 length vector
	%	with its covariance
	if (nargout > 2)
		[~,T] = normalize_data_isotropically(data_points);
		E = diag([1,2^-1,1,2^-1,2^-1,1]);
		% permutation matrix for interchanging 3rd and 4th
		% entries of a length-6 vector
		P34 = kron(diag([0,1,0]), [0 1;1 0]) + kron(diag([1,0,1]), [1 0; 0 1]);
		% 9 x 6 duplication matrix
		D3 = [1 0 0 0 0 0;
			    0 1 0 0 0 0;
			    0 0 1 0 0 0;
			    0 1 0 0 0 0;
			    0 0 0 1 0 0
			    0 0 0 0 1 0;
			    0 0 1 0 0 0;
			    0 0 0 0 1 0;
			    0 0 0 0 0 1];
			
		theta_guaranteed_covweightedNormalised = E \ P34*pinv(D3)*...
			inv(kron(T,T))'*D3*P34*E*theta_fastguaranteed;
		[S, thetaCovarianceMatrixNormalisedSpace] = compute_covariance_of_sampson_estimate(...
                     theta_fastguaranteed, data_points,covList);
		T_param2norm = T;
		ellipse_norm = theta_guaranteed_covweightedNormalised;
		cov_ell_norm = thetaCovarianceMatrixNormalisedSpace; 
	end
	
	% Use New method for computing the ellipse and its covariance
	geometricEllipseParameters = fromAlgebraicToGeometricParameters(theta_fastguaranteed);
	% Ellipse is computed under the form : a, b, x0, y0, theta
	rot = zeros(5,5);
	rot(1,3)=1; rot(2,4)=1; rot(3,1)=1; rot(4,2)=1; rot(5,5)=1;
	ellipse_param = transpose(rot*geometricEllipseParameters);
	if abs(imag(ellipse_param(5))) > 0
		geometricEllipseParameters
		theta_fastguaranteed
		save('data_points_extract.mat','data_points');
	end
	ellipse_param(5) = fix_angle(ellipse_param(5));
	% Sort the ellipse parameters so that if ellipse_pa
	if nargin > 1
		geoCov =  compute_covariance_of_geometric_parameters(theta_fastguaranteed, data_points,covList);

	else
		geoCov =  compute_covariance_of_geometric_parameters(theta_fastguaranteed, data_points);
	end
	cov_ell_param = rot*geoCov*rot';
	ellipse_param(1:2) = ellipse_param(1:2);
end


% This function adapt the angle so that it stands between -pi/2 and pi/2
function [angle_out] = fix_angle(angle_in);
if abs(imag(angle_in))>0
	angle_in = 0;
end
angle_out = mod(angle_in,2*pi);
if angle_out > pi
	angle_out = angle_out-pi;
end
end

