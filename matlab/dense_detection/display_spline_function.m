function [fig_handle] = display_spline_function(r_opt,r_proj,I_vect,I_pt,fig_name)
	nb_pts = 100;
	NB_R = length(r_opt);
	R_vect = tril(ones(NB_R));
        inv_R_vect = inv(R_vect);
	a_opt = inv_R_vect*r_opt;
	xx = r_opt(1):(r_opt(end)-r_opt(1))/nb_pts:r_opt(end);
	% Check if we can apply the monotone algorithm
	if a_opt(2:end)>0 
		pp = monotone_spline(r_opt,transpose(I_vect));
		yy = ppval(pp,xx);
	else
		yy = spline(r_opt,I_vect,xx);
	end
	I_proj = ppval(pp,r_proj);
	fig_handle = figure('Name',fig_name);
	plot(r_proj(:),I_pt(:),'+r');
	hold on;
	plot(r_opt,I_vect,'+g','MarkerSize',10,'LineWidth',3);
	plot(xx,yy,'-b','LineWidth',3);
	xlim([min(r_proj(:)),max(r_proj(:))]);
	title(fig_name);
	% Order the value to calculate a sliding window for standard deviation
	[I_pt_sorted,sorting_index] = sort(I_pt);
	r_proj_sorted = r_proj(sorting_index);
	%figure('Name','Sorted by intensities');
	%plot(I_pt_sorted,r_proj_sorted,'+k');
	%NB_DIV = 100;
	%nb_pt = length(I_pt_sorted);
	%r_proj_matrix = reshape(r_proj_sorted(1:(end-mod(nb_pt,NB_DIV))),floor(nb_pt/NB_DIV),NB_DIV);
	%figure('Name','Variance for each bin');
	%plot(var(r_proj_matrix));
	
end
