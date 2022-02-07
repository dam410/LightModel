function [data] = isocontour_detection(data,I,polys_2D,options)
	arguments
		data (1,1) struct
		I double
		polys_2D cell
		options.Mode (1,1) string = 'top-down'
		options.Setting (1,1) string = 'perspectivity'
		options.IsocontourDiscretisation (1,1) string = 'virtual'
		options.Display (1,1) string = 'off'
		options.N_circ (1,1) double = 1
		options.NB_R (1,1) double = 5
	end
	switch options.Mode
		case 'bottom-up'
			[pts,levels,curve_detected] = detect_isocontour_real(I,polys_2D,options.N_circ,options.Display);
			isocontour_curve_params = {};
			isocontour_pts = {};
			for i=1:size(polys_2D,2)
				isocontour_pts{i} = pts(i,:);
				isocontour_curve_params{i} = curve_detected(i,:);
			end
			isocontour.Points = isocontour_pts;
			isocontour.CurveParameters = isocontour_curve_params;
			data.isocontour = isocontour;
		case 'top-down'
                        [homog_p] = dense_isocontours_detection(I,data.K,polys_2D,...
                                'Mode',options.Setting,...
                                'Display',options.Display,...
                                'Nb_R',options.NB_R);
                        data = add_isocontours_to_homog(data,polys_2D,homog_p,options.N_circ);
		case 'top-down-circular'
			[homog_p] = dense_isocontours_detection(I,data.K,polys_2D,...
				'Mode','rotation',...
				'Display',options.Display,...
				'Nb_R',options.NB_R);
			data = add_isocontours_to_homog(data,polys_2D,homog_p,options.N_circ);
	end
		
end
