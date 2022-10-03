function [position] = getRectROI(I);
	h = figure('Name','Selection of the rectangle');
	imshow(I/255);
	roi = drawrectangle();
	position = roi.Position
	close(h);
end
