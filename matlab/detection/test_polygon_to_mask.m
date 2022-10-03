poly_2D = [100,400,400,500,100,-100;-200,-100,100,300,300,100];
n = 600;
m = 400;
I = zeros(n,m);
figure('Name','Display the points and the image boundary');
plot([1,1,n,n,1],[1,m,m,1,1],'-r');
hold on;
plot(poly_2D(1,:),poly_2D(2,:),'+b');
res_mask_1 = polygon_to_mask(poly_2D,I,0);
res_mask_2 = polygon_to_mask(poly_2D,I,10);

figure('Name','Display the mask 1, the points and the image boundary');
imshow(res_mask_1);
hold on;
plot([1,1,n,n,1],[1,m,m,1,1],'-r');
plot(poly_2D(1,:),poly_2D(2,:),'+b');

figure('Name','Display the mask 2, the points and the image boundary');
imshow(res_mask_2);
hold on;
plot([1,1,n,n,1],[1,m,m,1,1],'-r');
plot(poly_2D(1,:),poly_2D(2,:),'+b');


BW = poly2mask(poly_2D(1,:),poly_2D(2,:),m,n);
BW = bwmorph(BW,'shrink',5);

figure('Name','With direct matlab function');
imshow(BW);
hold on;
plot(poly_2D(1,:),poly_2D(2,:),'+b');

figure('Name','With only the contour points');
BW = bwmorph(BW,'remove');
imshow(BW);
hold on;
plot(poly_2D(1,:),poly_2D(2,:),'+b');





