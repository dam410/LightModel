# To use execute those lines in blender
import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.gen_file_scene(bpy,'test_X.json','test_X.png','/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic');
get_data_script.gen_file_scene(bpy,'test_noise.json','test_noise.png','/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic');
# Example
get_data_script.create_realdata_problem_replicate(bpy,'/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic_realdata_replicate');


get_data_script.generate_data_one_plane_auto(bpy,'/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic');

# Data for scene with two planes
importlib.reload(get_data_script);
get_data_script.generate_data_two_planes_auto(bpy,'/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic');




import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.place_plane(bpy,0.7,Vector([0,0,-1]))

import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.rotation_euler(get_data_script.rotation_matrix(1,0.4,-0.9))

import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.variatiate_angle(bpy,80,5,"/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic_exp_1");

import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.nb_plane_variation(bpy,20,"/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic_exp_1")

import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.distance_source_variation(bpy,30,10,"/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic_exp_2")

import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.noise_level_variation(bpy,6,10,"/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic_exp_2")


import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.angle_variation(bpy,80,5,"/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic_exp");

import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.distance_source_variation(bpy,40,10,"/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic_exp_2");



import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.place_plane(bpy,1.5,0,D.objects[1],1.0);


import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import create_scene_facets
importlib.reload(create_scene_facets)
scene = D.scenes[0];
cam = D.objects[0];
#create_scene_facets.create_random_grid('grid',8,5,scene,cam)
#create_scene_facets.create_random_grid('grid',3,3,scene,cam,0.00,0.8);
grid = create_scene_facets.create_random_grid('grid',3,3,scene,cam,10.00,0.8);


# To use execute those lines in blender
import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import get_data_script
importlib.reload(get_data_script);
get_data_script.gen_file_scene(bpy,'test_multifacet.json','test_multifacet.png','/home/dam/Documents/PostDoc_Damien/LightModel/data/synthetic');

import importlib
import sys
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
import create_scene_facets
importlib.reload(create_scene_facets)
scene = D.scenes[0];
cam = D.objects[0];


import importlib
import sys
sys.path.append('/home/mourad/Documents/PostDoc_Damien/Git/LightModel/blender/');
import script_endo_experiment
import create_scene_facets
import get_data_script
importlib.reload(script_endo_experiment);
importlib.reload(get_data_script);
importlib.reload(create_scene_facets);
script_endo_experiment.exp_source_camera_distance();


import importlib
import sys
sys.path.append('/home/mourad/Documents/PostDoc_Damien/Git/LightModel/blender/');
import script_endo_experiment
import create_scene_facets
importlib.reload(script_endo_experiment);
importlib.reload(create_scene_facets);
create_scene_facets.create_random_grid('test',2,2,bpy.data.scenes[0],bpy.data.objects[0],0.26,1.0);

import importlib
import sys
sys.path.append('/home/mourad/Documents/PostDoc_Damien/Git/LightModel/blender/');
import script_endo_experiment
import create_scene_facets
import get_data_script
importlib.reload(script_endo_experiment);
importlib.reload(get_data_script);
importlib.reload(create_scene_facets);
script_endo_experiment.exp_source_camera_one_scene();

