N = 30; 
pts_1 = rand(2,N);
pts_2 = rand(2,N);
figure('Name','Display segments ans their intersections');
plot([pts_1(1,:);pts_2(1,:)],[pts_1(2,:);pts_2(2,:)],'-b');
hold on;
for i = 1:(N-1)
	for j = i:N
		[pt] = seg_inter([pts_1(:,i),pts_2(:,i)],[pts_1(:,j),pts_2(:,j)]);
		if length(pt)>1
			plot(pt(1),pt(2),'+r');
		end
	end
end
