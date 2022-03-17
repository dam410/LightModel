function [N_out] = convention_normals(N,Xc,K,pt_vis)
	% Apply the convention for the normals
	pt_vis
	if ~isempty(pt_vis)
		P = inv(K)*[transpose(pt_vis);1];
		P = P/P(3);
		N_out = sign(transpose(P)*N).*N;
	else
		dot_NX = sum(Xc.*N(1:2,:),1) + N(3,:);
		N_out = sign(dot_NX).*N;
	end
end
