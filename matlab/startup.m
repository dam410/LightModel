setPath;
createHandleSymbolic;
mex -L"/usr/local/lib" -lapriltag -I"/usr/local/include/apriltag" mex/detect_apriltag.c
%addpath ~/Source/matlab_toolbox_cv/rvctools/common
%startup_rvc
