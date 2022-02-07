function [] = display_normals(all_planes)
	nb_plane = size(all_planes,2);
	Ns1 = cell2mat(all_planes(1,:));
	Ns2 = cell2mat(all_planes(2,:));
	figure('Name','Displaying all the estimated normal pairs');
	quiver3(zeros(1,nb_plane),zeros(1,nb_plane),zeros(1,nb_plane),...
		Ns1(1,:),Ns1(2,:),Ns1(3,:),...
		'MarkerMode','manual','Color','r','MaxHeadSize',0.2,...
		'Linewidth',2,'MarkerSize',5,'Autoscale','off','AutoScaleFactor',2);
	hold on;
	quiver3(zeros(1,nb_plane),zeros(1,nb_plane),zeros(1,nb_plane),...
		Ns2(1,:),Ns2(2,:),Ns2(3,:),...
		'MarkerMode','manual','Color','b','MaxHeadSize',0.2,...
		'Linewidth',2,'MarkerSize',5,'Autoscale','off','AutoScaleFactor',2);
end
