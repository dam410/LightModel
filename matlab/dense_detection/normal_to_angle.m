function [alpha,beta] = normal_to_angle(N)
        N = N/norm(N);
        beta = asin(N(3));
        if abs(cos(beta))<1e-10
                alpha = 0;
        else
                alpha = atan2(N(2)/cos(beta),N(1)/cos(beta));
        end
end
