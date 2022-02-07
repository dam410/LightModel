function [L_cell] = lines_2p_to_plucker(pts_A,pts_B)
	L_cell = cell(1,size(pts_A,2));
	for i=1:length(L_cell)
		L = [pts_A(:,i);1]*transpose([pts_B(:,i);1])-...
			[pts_B(:,i);1]*transpose([pts_A(:,i);1]);
		L_cell{i} = L;
	end
end
