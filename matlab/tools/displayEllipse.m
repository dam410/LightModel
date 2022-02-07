% Color can be a letter corresponding to the color or an 3 element-array representing RGB 
function [h] = displayEllipse(C,color,namedFigure)

if nargin<2
    color = 'g';
end
if nargin>=3
    namedFigure; hold on;
end

if (size(C,1) == 3) && (size(C,2) == 3)
    param = ellipse2param(C);
elseif (size(C,1) == 1) && (size(C,2) == 5)
    param = C;
else
    error('First argument in displayEllipse must be of size 1x5 or 3x3');
end
pts = ellipsepoints(param,100000);
if param(3)*param(4) < 0
	% On peut chercher un changement de signe de la dérivée 2 ??
	diff_x_pts = abs(pts(1,1:end)-pts(1,[end,1:(end-1)]));
	diff_y_pts = abs(pts(2,1:end)-pts(2,[end,1:(end-1)]));
	ind_valid = find(diff_x_pts<10 & diff_y_pts<10);
	diff_valid = ind_valid-[ind_valid(1),ind_valid(1:(end-1))];
	cut_curve = [1,find(diff_valid>1),length(ind_valid)];
	for k=1:(length(cut_curve)-1)
		hold on;
		inter_val = ind_valid(cut_curve(k):(cut_curve(k+1)-1));
		if length(color) == 1
			h = plot(pts(1,inter_val),pts(2,inter_val),color);
			set(h,'linewidth',2);
		else
			h = plot(pts(1,inter_val),pts(2,inter_val),'Color',color);
			set(h,'linewidth',2);
		end
	end
else 
	if length(color) == 1
		h = plot(pts(1,:),pts(2,:),color);
	else
		h = plot(pts(1,:),pts(2,:),'Color',color);
	end
end
set(h,'linewidth',2);
if nargin>=3
	hold off;
end
