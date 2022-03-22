# Script to generate images with different baseline and random orientation of the 3x3 tiles grid
# We suppose that scene has already been loaded with at least a light source called 'Lamp'
#	and a camera call 'Camera'
import math
from mathutils import *; from math import *; from numpy import *
import json
import itertools
import create_scene_facets
import get_data_script
import bpy

# Code for generating N random grid and acquisitions
def generate_samples(name_gen,output_folder,output_files_list,nb_sample,ang_deg,dis_avg):
	# Load scene
	scene = bpy.data.scenes[0];
	cam = bpy.data.objects[0];
	# Degree to radian conversion
	ang_var = math.pi*ang_deg/180.0;
	# For N images
	for i_sample in  range(0,nb_sample):
		# Generate the mesh
		create_scene_facets.create_random_grid('temp_mesh',2,3,scene,cam,ang_var,dis_avg);
		# Capture image
		i_str = "_"+"{:0>2d}".format(int(i_sample));
		generated_filename = name_gen + i_str;
		get_data_script.gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',output_folder);
		output_files_list.append(generated_filename+".json");
		# Remove the mesh
		bpy.ops.object.select_all(action='DESELECT');
		if bpy.app.version > (2,80,0):
			bpy.data.objects['temp_mesh'].select_set(True);
		else:
			bpy.data.objects['temp_mesh'].select = True;
		bpy.ops.object.delete();


def gen_sample(scene,cam,ang_var,name_gen,output_folder,output_files_list,dis_avg):
	# Generate the mesh
	create_scene_facets.create_random_grid('temp_mesh',2,3,scene,cam,ang_var,dis_avg);
	# Capture image
	i_str = "_"+"{:0>2d}".format(int(i_sample));
	generated_filename = name_gen + i_str;
	get_data_script.gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',output_folder);
	output_files_list.append(generated_filename+".json");
	# Remove the mesh
	bpy.ops.object.select_all(action='DESELECT');
	if bpy.app.version > (2,80,0):
		bpy.data.objects['temp_mesh'].select_set(True);
	else:
		bpy.data.objects['temp_mesh'].select = True;
	bpy.ops.object.delete();


# Generate experiment data
def exp_source_camera_distance():
	nb_dist = 20;
	bpy.data.objects.get('Lamp').location;
	output_folder = '/home/mourad/Documents/PostDoc_Damien/Git/LightModel/data/exp_synth_source_camera'
	exp_name = 'exp_synth_source_cam_dis_';
	out_list = list();
	scene = bpy.data.scenes[0];
	cam = bpy.data.objects[0];
	# Degree to radian conversion
	ang_var = math.pi*ang_deg/180.0;
	# Distance configuration with discretisation
	for i_d in range(0,nb_dist):
		bpy.data.objects.get('Lamp').location = Vector([i_d*0.5/nb_dist,0.0,0.0]);
		name_gen = exp_name + str(i_d)
		generate_samples(name_gen,output_folder,out_list,10,25,1.0);
	# Write a file with all the name
	f = open(output_folder+"/filenames_source_cam.txt", "a");
	for filename in out_list:
		f.write(filename + "\n");
	f.close();

def exp_source_camera_one_scene():
	nb_dist = 20;
	dis_avg = 1.0;
	ang_var = math.pi*45/180.0;
	bpy.data.objects.get('Lamp').location;
	output_folder = '/home/mourad/Documents/PostDoc_Damien/Git/LightModel/data/exp_synth_source_camera'
	exp_name = 'exp_synth_source_cam_one_scene_dis_';
	out_list = list();
	scene = bpy.data.scenes[0];
	cam = bpy.data.objects[0];
	# Degree to radian conversion
	# Generate one scene
	create_scene_facets.create_random_grid('temp_mesh',2,3,scene,cam,ang_var,dis_avg);
	# Capture image
	for i_d in range(0,nb_dist):
		bpy.data.objects.get('Lamp').location = Vector([i_d*0.5/nb_dist,0.0,0.0]);
		name_gen = exp_name + str(i_d)
		generated_filename = name_gen;
		get_data_script.gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',output_folder);
		out_list.append(generated_filename+".json");
	f = open(output_folder+"/filenames_source_cam_one_scene.txt", "a");
	for filename in out_list:
		f.write(filename + "\n");
	f.close();


