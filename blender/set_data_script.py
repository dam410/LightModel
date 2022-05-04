from mathutils import *;
from math import *;
from numpy import *
import json
import itertools
import math
import bpy

def read_scene(bpy,json_file,folder):
	f = open(folder+'/'+json_file,'r');
	scene_dict = json.loads(f.read());
	f.close();
	return scene_dict

def set_lamp_dict(bpy,lamp_dict):
	# Check and find the lamp if it exists otherwise create a lamp
	if bpy.data.objects.find('Lamp'):
		lamp = bpy.data.objects.get('Lamp');
		lamp.data.type = lamp_dict['type'].upper();
	else:
		lamp = bpy.data.lights.new(name='Lamp', type=lamp_dict['type'].upper());
	if lamp_dict['type']=='point':
		lamp.data.falloff_type = lamp_dict['data']['falloff_type'];
		lamp.data.quadratic_attenuation = lamp_dict['data']['quadratic_attenuation'];
	elif lamp_dict['type']=='spot':
		lamp.data.falloff_type = lamp_dict['data']['falloff_type'];
		lamp.data.quadratic_attenuation = lamp_dict['data']['quadratic_attenuation'];
		lamp.data.spot_size = lamp_dict['data']['spot_size'];
		lamp.data.spot_blend = lamp_dict['data']['spot_blend'];
	else:
		print('Import Error: Type of lamp not recognized');

def set_cam_dict(bpy,cam_dict):
	if bpy.data.objects.find('Camera'):
		cam = bpy.data.objects.get('Camera'):
	cam = bpy.data.cameras.new(name='Camera');
	cam.data.lens = cam_dict['focal'];
	cam.data.type = cam_dict['type'];
	cam.data.sensor_fit = cam_dict['sensor_fit'];
	cam.data.sensor_height = cam_dict['sx'];
	cam.data.sensor_width = cam_dict['sy'];
	cam.data.shift_x cam.data.shift_xcam_dict['x0'];
	cam.data.shift_y = cam_dict['y0'];

def remove_mesh():
        # Remove all the previous meshes
        bpy.ops.object.select_all(action='DESELECT');
        for obj in bpy.context.scene.objects:
                if obj.type in {'CURVE','SURFACE','MESHES','POINTCLOUD','VOLUME'}:
                        if bpy.app.version > (2,80,0):
                                obj.select_set(True);
                        else:
                                obj.select = True;
        bpy.ops.object.delete();

def set_mesh_dict(mesh_dict):
	mesh = bpy.data.meshes.new('imported_mesh');
	# Add the geometry to the mesh according to the data
	all_vert = mesh_dict['data']['vertices'];
	all_faces = list(tuple(x['vertices_index']) for x in mesh_dict['data']['polygons']);
	mesh.from_pydata(all_vert, [], all_faces);
	# Create the material
	for i_mat in range(0,len(mesh_dict['data']['materials'])):
		mat_dict = mesh_dict['data']['materials'][i_mat];
		mat = bpy.data.materials.new(name=('Material_'+str(i_mat)));
		if bpy.app.version > (2,80,0):
			dif_value = mat_dict['diffuse_intensity'];
			bpy.data.materials[i_mat].diffuse_color = Color((dif_value,dif_value,dif_value));
		else:
			bpy.data.materials[i_mat].diffuse_color = dif_value;
			bpy.data.materials[i_mat].specular_shader = mat_dict['specular_shader'];
			bpy.data.materials[i_mat].diffuse_shader = mat_dict['diffuse_shader'];
			bpy.data.materials[i_mat].ambient = mat_dict['ambient'];
		bpy.data.materials[i_mat].specular_intensity = mat_dict['specular_intensity'];
	for i_pol in range(0,mesh_dict['data']['polygons']):
		poly_dict = mesh_dict['polygons'][i_pol];
		poly = mesh.data.polygons[i_poly]
		poly.material_index = poly_dict['material_index'];

def set_pose_obj(obj_dict):
			
		
		
