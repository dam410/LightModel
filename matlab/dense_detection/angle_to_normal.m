function [N,J_N] = angle_to_normal(alpha,beta)
        N = [cos(alpha)*cos(beta);sin(alpha)*cos(beta);sin(beta)];
        J_N = [-cos(beta)*sin(alpha), -cos(alpha)*sin(beta);...
                cos(alpha)*cos(beta), -sin(alpha)*sin(beta);...
                     0, cos(beta)];
end