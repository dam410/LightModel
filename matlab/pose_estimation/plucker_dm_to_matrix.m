function [L] = plucker_dm_to_matrix(d,m)
	%L = [0,-d(1),-d(2),-d(3);d(1),0,-m(1),-m(2);d(2),m(1),0,-m(3);d(3),m(2),m(3),0];
	L = [[cross_antisym(m),d];[-transpose(d),0]];
end
