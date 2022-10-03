function [] = img_fusion(list_img,img_out)
	
	I = double(imread(list_img{1}));
	for i=2:length(list_img)
		I = I + double(imread(list_img{i}));
	end
	I = I/length(list_img);
	I = uint8(I);
	imwrite(I,img_out);
end
