function [rho] = camera_response(L)
	k = 10;
	rho = 255./(1+exp(-k*(2*L-1)));
end

