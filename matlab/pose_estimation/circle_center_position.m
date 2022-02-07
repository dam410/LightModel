function [X] = circle_center_position(K,Xc,N,d)
	X_ray = inv(K)*[Xc;1];
        lambda = -d/(transpose(N)*X_ray);
        X = lambda*X_ray;
end
