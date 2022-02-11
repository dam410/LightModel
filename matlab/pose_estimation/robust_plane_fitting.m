% Function to estimate common direction of planes intersecting at origin
function [L,err] = robust_plane_fitting(planes)
	Ns = cell2mat(planes);	
	M = Ns(1:3,:);
        A = M*transpose(M);
	b = robustfit(M,zeros(size(M,2),1));
        try
                [U,S,V] = svd(A);
        catch e
                disp('Problem with intersection of planes');
                V = eye(3);
                S = eye(3);
        end
        L = plucker_dm_to_matrix(V(:,3),zeros(3,1));
        err = S(3,3);
        G = zeros(size(M));
        figure;
        plot3([G(1,:);M(1,:)],[G(2,:);M(2,:)],[G(3,:);M(3,:)],'-b');
	hold on;
	plot3([0,b(1)],[0,b(2)],[0,b(3)],'-r');
	plot3([0,L(1)],[0,L(2)],[0,L(3)],'-g');
        axis equal;

	
end
