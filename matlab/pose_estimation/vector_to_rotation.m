function [H] = vector_to_rotation(h_vect,invK)
	N = angle_to_normal(h_vect(1),h_vect(2));
	base_R = null(N*transpose(N));
	if dot(cross(base_R(:,1),base_R(:,2)),N)<0
                R = [base_R(:,[2,1]),N];
        else
                R = [base_R(:,[1,2]),N];
        end
        %else
        %       R = [base_R(:,[1,2]),N];
        %end
        H = inv(R)*invK;
end
