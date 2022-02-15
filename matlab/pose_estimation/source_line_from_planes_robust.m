% Function to estimate common direction of planes intersecting at origin
function [L_mat,err] = source_line_from_planes_robust(planes)
	Ns = cell2mat(planes);	
	M = Ns(1:3,:);
	fitFcn = @(x) fitting_source_line(x);
	distFcn = @(x,y) err_L1(x,y);
	[L,inliersIdx] = ransac(transpose(M),fitFcn,distFcn,2,1e-1);
        [L_mat] = plucker_dm_to_matrix(transpose(L),zeros(3,1));
        err = sum(err_L1(L,transpose(M)));
        %G = zeros(size(M));
        %figure;
        %plot3([G(1,:);M(1,:)],[G(2,:);M(2,:)],[G(3,:);M(3,:)],'-r');
	%hold on;
        %plot3([G(1,inliersIdx);M(1,inliersIdx)],...
	%	[G(2,inliersIdx);M(2,inliersIdx)],...
	%	[G(3,inliersIdx);M(3,inliersIdx)],'-b');
	%plot3([0,L(1)],[0,L(2)],[0,L(3)],'-g');
        %axis equal;
end

function [L] = fitting_source_line(D)
	A = transpose(D)*D;
        try
                [U,S,V] = svd(A);
        catch e
                disp('Problem with intersection of planes');
                V = eye(3);
                S = eye(3);
        end
	L = transpose(V(:,3));
end

function [res] = err_L1(L,D)
	res = abs(D*transpose(L));
end
