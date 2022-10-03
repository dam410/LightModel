I = zeros(800,1400);
size_I = size(I);
n = size_I(1);
m = size_I(2);
% Two cases intersection or no intersection
sx = floor(size_I(1)/2);
sy = floor(size_I(2)/2);
pt_from = [sx;sy];
pt_ext = [3*sx;4*sy];
pt_int = floor([1.5*sx;1.5*sy]);
pt_ext_2 = floor([2.5*sx;0.5*sy]);
[inboundcheck,pt_inter] = in_bc(pt_int,size_I,pt_from)
[inboundcheck,pt_inter] = in_bc(pt_ext,size_I,pt_from)
[inboundcheck,pt_inter_2] = in_bc(pt_ext_2,size_I,pt_from)
figure('Name','Check my boundary function');
plot([1,1,n,n,1],[1,m,m,1,1],'-r');
hold on;
plot([pt_from(1),pt_ext(1)],[pt_from(2),pt_ext(2)],'-b');
plot([pt_from(1),pt_int(1)],[pt_from(2),pt_int(2)],'-b');
plot([pt_from(1),pt_ext_2(1)],[pt_from(2),pt_ext_2(2)],'-b');
plot(pt_inter(1),pt_inter(2),'+g');
plot(pt_inter_2(1),pt_inter_2(2),'+g');

