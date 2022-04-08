function [h_vect] = homography_to_vector_affine(H)
        H_aff = H(1:2,:);
        h_vect = H_aff(:);
end
