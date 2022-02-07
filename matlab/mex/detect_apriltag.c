#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <apriltag.h>
#include <tag16h5.h>
#include <apriltag_pose.h>
#include <math.h>
#include "matrix.h"
#include "common/getopt.h"
#include "common/image_u8.h"
#include "common/image_u8x4.h"
#include "common/pjpeg.h"
#include "common/zarray.h"

/** [cell_R,cell_t,cell_id] = detect_apriltag(img_filepath_pnm,K,tagSize)
 * Input Parameters :
 *  img_filepath_pnm : filepath to image in pnm without distorsion
 *  K : Intrinsic Calibration Matrix 
 *  tagSize : size of the tag in meter
 * Output Parameters :
 *  cell_R : cell array containing the orientation matrix (3x3) of the detected markers
 *  cell_t : cell array containing the translation vector (3x1) of the detected markers
 *  array_id : Array containing the id of the markers (int)
 *  array_err : Array containing the error 
*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Copy input variables */
	double *K = (double*) mxGetPr(prhs[1]);
	double *tagSize = (double*) mxGetPr(prhs[2]);
	/* Get the string of the filepath */
	char *input_buf;
	size_t buflen;
	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
	input_buf = mxArrayToString(prhs[0]);	

	/* Prepare output matrices */
	mxArray *cell_R;
	mxArray *cell_t;
	mxArray *array_id;
	mxArray *array_err;


	/* Detection of the markers */	
        image_u8_t* im = image_u8_create_from_pnm(input_buf);
        apriltag_detector_t *td = apriltag_detector_create();
        apriltag_family_t *tf = tag16h5_create();
        apriltag_detector_add_family_bits(td, tf,0);
        zarray_t *detections = apriltag_detector_detect(td, im);

	/* Create temporary array */

	int nb_detection = zarray_size(detections);

	// Set output
	//cell_R = plhs[0] = mxCreateCellMatrix(1, nb_detection);
	//cell_t = plhs[1] = mxCreateCellMatrix(1, nb_detection);
	//array_id = plhs[2] = mxCreateDoubleMatrix(1, nb_detection, 0);
	//array_err = plhs[3] = mxCreateDoubleMatrix(1 ,nb_detection, 0);
	cell_R = plhs[0] = mxCreateCellMatrix(1, nb_detection);
	cell_t = plhs[1] = mxCreateCellMatrix(1, nb_detection);
	array_id = plhs[2] = mxCreateDoubleMatrix(1, nb_detection, 0);
	array_err = plhs[3] = mxCreateDoubleMatrix(1 ,nb_detection, 0);


	double *id_,*err_;
	id_ = mxGetPr(array_id);
	err_ = mxGetPr(array_err);

        for (int i = 0; i < nb_detection; i++) {
                apriltag_detection_t *det;
                zarray_get(detections, i, &det);
                apriltag_detection_info_t info;
                info.det = det;
                info.tagsize = tagSize[0];
                info.fx = K[0];
                info.fy = K[4];
                info.cx = K[6];
                info.cy = K[7];
                apriltag_pose_t pose;
                double err = estimate_tag_pose(&info, &pose);
		//printf("Detection number=%d with error=%f\n",det->id,err);
		id_[i] = det->id;
		err_[i] = err;

		//mxArray *R_ = mxCreateDoubleMatrix(3, 3, 0);
		//mxArray *t_ = mxCreateDoubleMatrix(3, 1, 0);

		//mxSetCell(cell_t, i, t_);
		//mxSetCell(cell_R, i, R_);

		//memcpy((double*) mxGetPr(R_),(double*) R->data,sizeof(double)*9);
		//memcpy((double*) mxGetPr(t_),(double*) t->data,sizeof(double)*3);

		// Remove previous allocated data
		//mxDestroyArray(mxGetCell(cell_R,i));
		//mxDestroyArray(mxGetCell(cell_t,i));
		mxSetCell(cell_t, i, mxCreateDoubleMatrix(3, 1, 0));
		mxSetCell(cell_R, i, mxCreateDoubleMatrix(3, 3, 0));

		memcpy((double*) mxGetPr(mxGetCell(cell_R,i)),(double*) pose.R->data,sizeof(double)*9);
		memcpy((double*) mxGetPr(mxGetCell(cell_t,i)),(double*) pose.t->data,sizeof(double)*3);

        }
	apriltag_detections_destroy(detections);

        //printf("After getting all pose data\n");


	
        /* Cleanup and memory freeing */
	apriltag_detector_clear_families(td);
        apriltag_detector_destroy(td);
        tag16h5_destroy(tf);
	image_u8_destroy(im);
	mxFree(input_buf);

        //printf("After freeing memory\n");
}
