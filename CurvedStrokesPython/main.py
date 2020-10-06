import json
import tools
import sys
import os
import numpy as np
import pre_processing
import sketch_dependency_graph
import sketch_wires
from camera import Camera
#import polyscope as ps
from copy import deepcopy

MAX_NB_SKETCHES = 10

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-file_name", default="toto.json", help="")
args = parser.parse_args()
file_name = args.file_name

if file_name == "toto.json":
	import inspect
	current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
	parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
	file_name = os.path.join(parent_dir, "StraightStrokesMatLab", "folder_save",
							 "student8", "house", "view1", "student8_house_bestScore_full.json")
print("input_file:", file_name)

# read full_JSON
fh = open(file_name, "r")
file_content = json.load(fh)
fh.close()
designer_name = file_content["designer"]
object_name = file_content["object_name"]
strokes_topology = file_content["strokes_topology"]
output_strokes_topology = deepcopy(file_content["strokes_topology"])
output_strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology, use_rdp=False)
intersections = file_content["intersections"]
cam_params = file_content["cam_param"]
cam = Camera(proj_mat=file_content["cam_param"]["P"],
			 focal_dist=file_content["cam_param"]["f"],
			 fov=file_content["cam_param"]["fov"],
			 t=file_content["cam_param"]["t"],
			 view_dir=file_content["cam_param"]["view_dir"],
			 principal_point=file_content["cam_param"]["principal_point"],
			 rot_mat=file_content["cam_param"]["R"],
			 K=file_content["cam_param"]["K"],
			 cam_pos=file_content["cam_param"]["C"],
			 vanishing_points_coords=file_content["cam_param"]["vp_coord"])
cam.compute_inverse_matrices()

clusters, bbox, old_new_stroke_mapping = pre_processing.pre_process(file_name, designer_name, object_name,
																	strokes_topology, intersections, cam_params)
bbox_diag = np.linalg.norm(bbox[:3] - bbox[3:])
tools.scale_bbox(bbox, 0.5)

dep_graph_template = sketch_dependency_graph.init_dep_graph(strokes_topology)
#ps.init()
dep_graph_template.median_length_3d = np.median(list(dep_graph_template.stroke_lengths_3d.values()))
# for clipping, get median length in 2D
stroke_lengths = [strokes_topology[dep_node.stroke_id]["linestring"].length
				  for dep_node in dep_graph_template.dependency_nodes]
median_length_2D = np.median(stroke_lengths)

for cluster in clusters:
	dep_graph = deepcopy(dep_graph_template)
	global_inter_id = [0]
	sketch_dependency_graph.append_to_dep_graph(dep_graph, cluster, strokes_topology)
	extreme_intersections_distances_per_stroke = []
	total_length_2d = 0.0
	for dep_node_id, dep_node in enumerate(dep_graph.dependency_nodes):
		dep_graph.dependency_nodes[dep_node_id].length_2d = strokes_topology[dep_node.stroke_id]["linestring"].length
		total_length_2d += strokes_topology[dep_node.stroke_id]["linestring"].length
		if dep_node.is_assigned:
			extreme_intersections_distances_per_stroke.append([])
			continue
		s_id = dep_node.stroke_id
		inter_set = strokes_topology[s_id]["curve_intersections"]
		arc_params = [np.mean(inter["curved_t_param"]) for inter in inter_set]
		extreme_intersections_distances_per_stroke.append(
			np.max(arc_params) - np.min(arc_params))
	dep_graph.total_length_2d = total_length_2d
	#print("cluster: ", cluster)
	for s_id in cluster:
		dep_node_id = dep_graph.stroke_id_to_dep_node_id[s_id]
		dep_node = dep_graph.dependency_nodes[dep_node_id]
		if dep_node.is_assigned:
			continue
		#print("s_id: ", s_id)
		#print("dep_node_id: ", dep_node_id)
		# get contexts
		contexts_3d = dep_graph.get_contexts_3d(s_id)
		if len(contexts_3d) == 0:
			print("Warning: no 3D context")
		all_context_curves = []
		for contexts_3d_id, context_3d in enumerate(contexts_3d):
			if len(context_3d.sketch_versions) == 0:
				print("context_3d.sketch_versions empty")
				continue
			# get 3D intersections
			intersections_3d = sketch_wires.get_intersections(s_id, context_3d,
															  strokes_topology,
															  dep_graph, cam,
															  global_inter_id)
			# get candidate curves
			candidate_curves, intersection_sets = sketch_wires.get_candidate_curves(
				s_id, intersections_3d, strokes_topology, cam, bbox_diag)

			# remove overly foreshortened curves
			candidate_curves, intersection_sets = sketch_wires. \
				filter_out_foreshortened_lines(candidate_curves,
											   intersection_sets,
											   dep_graph,
											   s_id,
											   strokes_topology,
											   median_length_2D,
											   bbox)

			for cand_curve_id, cand_curve in enumerate(candidate_curves):
				dep_inter_set = [intersections_3d[inter_id] for inter_id in
								 intersection_sets[cand_curve_id]]
				all_context_curves.append(
					[cand_curve, dep_inter_set, context_3d.sketch_versions])
			sketch_wires.project_tangents_to_3d(all_context_curves, dep_graph,
												s_id, cam)

		# merge over contexts
		merged_context_curves = sketch_wires.merge_over_contexts_straight_lines(all_context_curves, bbox_diag, s_id)
		# update_dep_graph
		dep_graph.update_dep_graph(merged_context_curves, dep_node_id)
		# update_stroke_score_function
		dep_graph.update_strokes_score_function([dep_node_id], strokes_topology,
												extreme_intersections_distances_per_stroke,
												cam)
		# PRUNING
		# best_candidate pruning
		if not strokes_topology[s_id]["is_ellipse"] and \
				len(dep_graph.dependency_nodes[dep_node_id].candidate_nodes.keys()) > 5:
			delete_sketch_version = dep_graph.accept_best_cand_line(dep_node_id)
			dep_graph.prune_sketch_versions(delete_sketch_version)

		# update scores for a set of strokes
		# score_strokes are the dep_node_ids !!!
		score_strokes = dep_graph.stroke_id_to_dep_node_id[
			dep_node.in_strokes].tolist()
		dep_graph.update_strokes_score_function(score_strokes, strokes_topology,
												extreme_intersections_distances_per_stroke,
												cam)
		# early assignment pruning
		score_strokes.append(dep_node_id)
		delete_sketch_version = dep_graph.early_assignment(score_strokes)
		dep_graph.prune_sketch_versions(delete_sketch_version)

		# max_nb sketch_versions pruning
		best_sketch_versions = dep_graph.get_n_best_sketch_versions(
			MAX_NB_SKETCHES)
		delete_sketch_version = []
		for sketch_version in dep_graph.sketch_versions_reference.keys():
			if not (sketch_version in best_sketch_versions):
				delete_sketch_version.append(sketch_version)
		dep_graph.prune_sketch_versions(delete_sketch_version)

		dep_graph.median_length_3d = np.median(list(dep_graph.stroke_lengths_3d.values()))
	#dep_graph.display_sketch_versions()
	#ps.show()
	best_curves, best_intersections = dep_graph.get_best_sketch()
	# backproject to output_strokes_topology
	for s_id in cluster:
		#print(s_id)
		dep_node_id = dep_graph.stroke_id_to_dep_node_id[s_id]
		curve = best_curves[dep_node_id]
		old_s_id = strokes_topology[s_id]["old_index"]
		#print("old_s_id: ", old_s_id)
		affected_stroke_ids = old_new_stroke_mapping[np.argwhere(old_new_stroke_mapping[:, 1] == old_s_id).flatten().tolist()][:, 0].tolist()
		affected_stroke_ids.append(old_s_id)
		#print(affected_stroke_ids)
		for aff_s_id in affected_stroke_ids:
			backprojected_curve = np.array(cam.lift_polyline_to_plane(output_strokes_array[aff_s_id],
																	  curve.plane_point,
																	  curve.plane_normal))
			#print("aff_s_id: ", aff_s_id)
			if len(backprojected_curve) == 0:
				#print("geometry null")
				continue
			if np.any(backprojected_curve < bbox[:3]) or np.any(backprojected_curve > bbox[3:]):
				#print("out of bbox")
				continue
			output_strokes_topology[aff_s_id]["points3D"] = backprojected_curve.tolist()
			output_strokes_topology[aff_s_id]["depth_assigned"] = True

# put curve information back to the json-file
fh = open(file_name, "r")
file_content = json.load(fh)
fh.close()
file_content["strokes_topology"] = output_strokes_topology
output_file_name = file_name
output_file_name = output_file_name.split("/")[-1]
output_file_name = output_file_name.split(".json")[0]
output_file_name += "_curves.json"
print("output_file: ", output_file_name)

if not os.path.exists("folder_save"):
	os.mkdir("folder_save")
fh = open(os.path.join("folder_save", output_file_name), "w")
json.dump(file_content, fh, indent=4)
fh.close()
