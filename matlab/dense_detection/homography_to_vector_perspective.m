function [h_vect] = homography_to_vector_perspective(H)
	H = H./H(3,3);
	h_vect = transpose(H(1:8));
end
