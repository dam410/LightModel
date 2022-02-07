For each case we have its number, and the availabe input parameter

case_X : X is the corresponding case
Xp : X is the number of scene plane in used, it can be 1, 2 or n

dssp : Distance of the Source to the Scene Plane
	* dssp is a cell, its length is the number of scene planes
	* each cell contains a cell whose the length is the number of solution
		* a solution is a scalar, distance between source and scene plane

psp : Pose of the Scene Plane
	* psp is a cell, its length is the number of scene planes
	* each cell contains a cell whose the length is the number of solution
		* a solution is a 4x1 vector reprensenting pose of a plane

ps : Position of the Source
	* ps is a cell, its length is the number of solution for the source position
	* Each cell contains a solution, it can be depending on the cases :
		* a 3x1 vector, the position of the source estimated
		* a 4x4 matrix, the plucker matrix of the estimated line
		passing through the source
		* a 4x1 vector, the normal and the distance to the origin of the
		estimated plane containing the source 

case_2_1p_dssp_psp means this is the second case, the position of the source
is the only unknown


