function [Nd] = estimate_plane(P);
	% Augment the vector to homogeneous coordinates
	P = [P;ones(1,length(P))];
	[U,S,V] = svd(P*transpose(P));
        Nd = V(:,4);
        % Normalize N and check that its sign is opposite 
	%	to the direction of the points
        Nd = Nd/norm(Nd(1:3))*sign(Nd(4));
end
