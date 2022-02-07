function [dssp] = calculate_dssp(ps,psp)
        nb_plane = length(psp);
        dssp = cell(1,nb_plane);
        for i_p = 1:nb_plane
                dssp{i_p} = {-(psp{i_p}{1}(4)+transpose(psp{i_p}{1}(1:3))*ps{1})};
        end
end

