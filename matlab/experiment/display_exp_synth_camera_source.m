load('results_exp_synth_source_camera.mat');

% Cleaning the results data 0 -> Inf
% This happens for ||S|| = 0
% 	and ||S|| = 17/20 at the 4th sample
err_orient_global_data = err_orient_global_all;
err_orient_global_data(isinf(err_orient_global_data)) = NaN;
err_orient_global_avg = mean(err_orient_global_data,2,'omitnan');

err_orient_colocalised_avg = mean(err_orient_colocalised_all,2);

err_s_global_avg = Inf*ones(1,20);
err_s_global_avg([2:16,18:20]) = mean(err_s_global_all([2:16,18:20],:),2);
err_s_global_avg(17) = mean(err_s_global_all(17,[1,2,3,5]));

% Display the results on a single curve
fig = figure('Name','Estimated Scene Plane Normal Error (degree)');
plot((0:19)/20,err_orient_global_avg,'-r',(0:19)/20,err_orient_colocalised_avg,'-b');
legend('Case H','Case G^*');
ylabel('Error on scene plane orientation (degree)');
xlabel('Distance between PLS and camera (meter)');
ax = gca;
ax.FontSize = 24;
fig.PaperPositionMode = 'manual';
orient(fig,'landscape');
saveas(fig,['~/Documents/PostDoc_Damien/LightModel/doc/results_exp_synth_camera_source.pdf']);


figure('Name','Estimated PLS Error (meter)');
plot((0:19)/20,err_s_global_avg,'-r');


