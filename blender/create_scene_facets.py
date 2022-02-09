import importlib
import sys
import numpy as np
import math
import bpy
sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
#import get_data_script

# In this file we will create a surface with various plane and orientation to study the curve

# Creating a subplane in a frustrum which image is in position (x0,y0,x1,y1)

# Getting the camera matrix in Matlab compatible format
def get_camera_matrices(cam,scene):
	focal = cam.data.lens;
	sx = cam.data.sensor_height;
	sy = cam.data.sensor_width;
	x0 = cam.data.shift_x;
	y0 = cam.data.shift_y;
	h = scene.render.resolution_y;
	l = scene.render.resolution_x;
	scale = scene.render.resolution_percentage / 100.0 * l / h;
	pixel_aspect_ratio = scene.render.pixel_aspect_x/scene.render.pixel_aspect_y;
	if cam.data.sensor_fit=='VERTICAL':
		su = 2.0*h/l*focal*scale/sx/pixel_aspect_ratio;
		sv = 2.0*h/l*focal*scale/sy;
	else:
		su = 2.0*h/l*focal*scale/sx;
		sv = 2.0*h/l*focal*scale/sy;
	H_im = np.matrix([[-(l-1)/2,0,l/2+0.5],[0,-(h-1)/2,h/2+0.5],[0,0,1]]);
	return H_im*np.matrix([[-sv, 0, y0 ], [0, su, x0],[ 0, 0, 1]]);

# Create a subplane at distance 1 of the 
def create_subplane(ray_00,ray_01,ray_11,ray_10,ray_center,d_cam,N):
	# Calculate the position in space of the central point at a fixed distance of the camera
	pt_center = d_cam*np.array(ray_center)/math.sqrt(np.dot(ray_center,ray_center));
	# Deduce the plane
	d_plane = -np.dot(pt_center,N);
	# Calculate the position of the points
	pt_00 = -d_plane*np.array(ray_00)/np.dot(ray_00,N);
	pt_10 = -d_plane*np.array(ray_10)/np.dot(ray_10,N);
	pt_01 = -d_plane*np.array(ray_01)/np.dot(ray_01,N);
	pt_11 = -d_plane*np.array(ray_11)/np.dot(ray_11,N);
	# To change into tuple list tuple(pt_01.tolist()[0])
	vertices = [tuple(pt_00),tuple(pt_01),tuple(pt_11),tuple(pt_10)];
	return vertices;

# Control the material of the object
def set_material(obj,k_d,k_s,color = (1.0,1.0,1.0)):
	mat = bpy.data.materials.new(name="MaterialName");
	mat.specular_intensity = k_s;
	if bpy.app.version[0] < 3:
		mat.diffuse_intensity = 1.0;
	mat.diffuse_color = (color[0]*k_d,color[1]*k_d,color[2]*k_d,0.0);
	obj.data.materials.append(mat);

# This function takes two interval in radian to generate a random normal vector
def random_normal(a_inter,b_inter):
	a = np.random.rand()*(a_inter[1]-a_inter[0])+a_inter[0];
	b = np.random.rand()*(b_inter[1]-b_inter[0])+b_inter[0];
	print([math.sin(a)*math.cos(b),math.sin(a)*math.sin(b),math.cos(a)]);
	return [math.sin(a)*math.cos(b),math.sin(a)*math.sin(b),math.cos(a)];

def random_distance(d_inter):
	return np.random.rand()*(d_inter[1]-d_inter[0])+d_inter[0];

# Create a meshgrid which subdivise the image and calculate the rays at the intersection
# then create small planes with the rays given a distance to the camera and a normal
def create_random_grid(name,nb_l,nb_h,scene,cam,ang_var,d_var):
	h = scene.render.resolution_y;	
	l = scene.render.resolution_x;
	K = get_camera_matrices(cam,scene);
	all_vert = [];
	all_faces = [];
	count = 0;
	for n_l in range(0,nb_l):
		for n_h in range(0,nb_h):
			# Calculate the 4 
			pt_00 = (np.matrix([n_l/nb_l*l,n_h/nb_h*h,1.0])*(np.linalg.inv(K)).transpose()).tolist()[0];
			pt_01 = (np.matrix([n_l/nb_l*l,(n_h+1)/nb_h*h,1.0])*(np.linalg.inv(K)).transpose()).tolist()[0];
			pt_11 = (np.matrix([(n_l+1.0)/nb_l*l,(n_h+1.0)/nb_h*h,1.0])*(np.linalg.inv(K)).transpose()).tolist()[0];
			pt_10 = (np.matrix([(n_l+1.0)/nb_l*l,n_h/nb_h*h,1])*(np.linalg.inv(K)).transpose()).tolist()[0];
			pt_center = (np.matrix([(n_l+0.5)/nb_l*l,(n_h+0.5)/nb_h*h,1])*(np.linalg.inv(K)).transpose()).tolist()[0];
			vertices = create_subplane(pt_00,pt_10,pt_11,pt_01,pt_center,-random_distance([1.2-d_var/2.0,1.2+d_var/2.0]),random_normal([-ang_var,ang_var],[-math.pi/2,math.pi/2]));
			all_vert.extend(vertices);
			all_faces.append(tuple([x for x in range(count,count+4)]))
			count = count + 4;
	mesh = bpy.data.meshes.new(name+'_mesh');	
	mesh.from_pydata(all_vert, [], all_faces);
	obj = bpy.data.objects.new(name, mesh);
	# As we placed our mesh relative to the camera, we apply to the object the same transformation of the camera
	obj.location = cam.location;
	obj.rotation_euler = cam.rotation_euler;
	# Set the material properties
	set_material(obj,0.8,0);
	# Attach the object ot the scene
	if bpy.app.version > (2,80,0):
		bpy.context.collection.objects.link(obj);
	else:
		scene.objects.link(obj);
	return obj;

#create_random_grid('grid',3,3,scene,cam)


