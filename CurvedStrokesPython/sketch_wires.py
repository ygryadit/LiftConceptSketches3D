import sketch_dependency_graph
import sys
import numpy as np
from more_itertools import distinct_combinations
import tools_3d
from shapely.geometry import Point
from copy import deepcopy
from math import acos
import time
import tools
import networkx as nx

def get_tangents_2d(intersections):
	tangents_curves = []
	tangents_straight_lines = []
	for inter in intersections:

		curve_inter_coords = inter["curved_linestring"].interpolate(0.5, normalized=True)
		inter_coords_t = inter["curved_linestring"].project(curve_inter_coords)
		for p_id, p in enumerate(inter["curved_linestring"].coords):
			p_t = inter["curved_linestring"].project(Point(p))
			if p_t > inter_coords_t:
				tangents_curves.append(np.array([inter["curved_linestring"].coords[max(p_id-1, 0)],
												 inter["curved_linestring"].coords[max(p_id, 1)]]))
				break
		straight_inter_coords = inter["straight_linestring"].interpolate(0.5, normalized=True)
		inter_coords_t = inter["straight_linestring"].project(straight_inter_coords)
		for p_id, p in enumerate(inter["straight_linestring"].coords):
			p_t = inter["straight_linestring"].project(Point(p))
			if p_t > inter_coords_t:
				tangents_straight_lines.append(np.array([inter["straight_linestring"].coords[max(p_id-1, 0)],
														 inter["straight_linestring"].coords[max(p_id, 1)]]))
				break
	return tangents_curves, tangents_straight_lines


def get_intersections(s_id, context_3d, strokes_topology, dep_graph, camera,
					  global_inter_id):
	intersections = []
	intersections_2d = strokes_topology[s_id]["curve_intersections"]
	tangents_curves, tangents_straight_lines = get_tangents_2d(intersections_2d)

	for prev_cand in zip(context_3d.dep_node_ids, context_3d.cand_node_ids):
		prev_dep_node = dep_graph.dependency_nodes[prev_cand[0]]

		line_3d = prev_dep_node.candidate_nodes[prev_cand[1]].geometry
		for inter_id, inter in enumerate(intersections_2d):
			if dep_graph.stroke_id_to_dep_node_id[inter["straight_line_idx"]] != prev_cand[0]:
				continue
			inter_coords = np.array(inter["curved_linestring"].interpolate(0.5, normalized=True))
			if not prev_dep_node.is_curve:
				vec_3d = (line_3d[-1] - line_3d[0])
				vec_3d /= np.linalg.norm(vec_3d)
				lifted_inter = camera.lift_point_close_to_line(
					inter_coords, line_3d[0], vec_3d)
			else:
				curve_plane_point = \
					prev_dep_node.candidate_nodes[prev_cand[1]].plane_point
				curve_plane_normal = \
					prev_dep_node.candidate_nodes[prev_cand[1]].plane_normal
				lifted_inter = camera.lift_point_to_plane(
					inter_coords, curve_plane_point, curve_plane_normal)
			dep_inter = sketch_dependency_graph.Intersection(
				inter_id=global_inter_id[0],
				coords_3d=lifted_inter,
				coords_2d=inter_coords,
				stroke_ids=[inter["straight_line_idx"], inter["curved_line_idx"]],
				tangents_2d=[tangents_straight_lines[inter_id], tangents_curves[inter_id]],
				inter_params=[[],[]],
				mid_inter_param=[np.mean(inter["straight_t_param"]), np.mean(inter["curved_t_param"])],
				acc_radius=inter["accuracy_radius"],
				adjacent_inter_ids=[])
			global_inter_id[0] += 1
			intersections.append(dep_inter)
	return intersections

def get_candidate_curves(s_id, intersections_3d, strokes_topology, camera, bbox_diag):

	candidate_curves = []
	# for each candidate line, include at least the intersections which were used
	# to create the candidate_line
	creation_intersections = []
	intersection_normals = [[] for i in range(len(intersections_3d))]

	#print("s_id: ", s_id)
	#start_time = time.clock()
	for inter_id, inter in enumerate(intersections_3d):
		if inter.coords_3d is None:
			continue
		# include the three major planes
		for i in range(3):
			cand_curve = tools_3d.Curve3D()
			plane_point = inter.coords_3d
			plane_normal = np.zeros(3)
			plane_normal[i] = 1.0
			#geometry = camera.lift_polyline_to_plane(strokes_topology[s_id]["stroke"],
			#										 plane_point, plane_normal)
			#print("geometry comparison")
			#print(np.array(geometry))
			intersection_normals[inter_id].append(plane_normal)
			geometry = camera.lift_polyline_to_plane_vectorized(strokes_topology[s_id]["stroke"],
																plane_point, plane_normal)
			#print(geometry)
			#sys.exit()
			#if len(geometry) == 0:
			#	continue
			cand_curve.geometry = geometry
			cand_curve.plane_point = plane_point
			cand_curve.plane_normal = plane_normal
			candidate_curves.append(cand_curve)
			creation_intersections.append([inter_id])

		if strokes_topology[s_id]["is_ellipse"]:
			continue
		curr_s_id = np.argwhere(np.array(inter.stroke_ids) == s_id).flatten()[0]
		other_s = strokes_topology[inter.stroke_ids[1-curr_s_id]]
		if other_s["primitive_type"] == 0:
			# add scaffold planes
			for plane in other_s["planes"]:
				plane_point = plane["plane_point"]
				plane_normal = plane["plane_normal"]
				if np.isclose(np.linalg.norm(plane_normal), 0.0):
					continue

				used_normals = np.array(intersection_normals[inter_id])
				if len(used_normals) > 0:
					used_dot = 1.0 - np.abs(np.dot(used_normals, plane_normal))
					if np.any(used_dot < np.deg2rad(0.1)/np.pi):
						continue
				intersection_normals[inter_id].append(plane_normal)

				geometry = camera.lift_polyline_to_plane_vectorized(strokes_topology[s_id]["stroke"],
														 plane_point, plane_normal)
				if len(geometry) == 0:
					continue
				cand_curve = tools_3d.Curve3D()
				cand_curve.geometry = geometry
				cand_curve.plane_point = plane_point
				cand_curve.plane_normal = plane_normal
				candidate_curves.append(cand_curve)
				creation_intersections.append([inter_id])

	#print("len(candidate_curves)")
	#print(len(candidate_curves))
	# finally, add planes formed by triplets of intersections
	for comb in distinct_combinations(range(len(intersections_3d)), 3):
		if strokes_topology[s_id]["is_ellipse"]:
			continue
		if intersections_3d[comb[0]].coords_3d is None or \
				intersections_3d[comb[1]].coords_3d is None or \
				intersections_3d[comb[2]].coords_3d is None:
			continue
		plane_point = np.array(intersections_3d[comb[0]].coords_3d)
		vec_1 = np.array(intersections_3d[comb[1]].coords_3d) - plane_point
		if np.isclose(np.linalg.norm(vec_1), 0.0):
			continue
		vec_1 /= np.linalg.norm(vec_1)
		vec_2 = np.array(intersections_3d[comb[2]].coords_3d) - plane_point
		if np.isclose(np.linalg.norm(vec_2), 0.0):
			continue
		vec_2 /= np.linalg.norm(vec_2)
		plane_normal = np.cross(vec_1, vec_2)
		if np.isclose(np.linalg.norm(plane_normal), 0.0):
			continue
		plane_normal /= np.linalg.norm(plane_normal)

		# check if similar normal already used by one of the 3 intersections
		normal_already_used = False
		for i in comb:
			used_normals = np.array(intersection_normals[i])
			if len(used_normals) > 0:
				used_dot = 1.0 - np.abs(np.dot(used_normals, plane_normal))
				if np.any(used_dot < np.deg2rad(0.1)/np.pi):
					normal_already_used = True
					break
		if normal_already_used:
			continue
		else:
			for i in comb:
				intersection_normals[i].append(plane_normal)

		geometry = camera.lift_polyline_to_plane_vectorized(strokes_topology[s_id]["stroke"],
												 plane_point, plane_normal)
		if len(geometry) == 0:
			continue
		cand_curve = tools_3d.Curve3D()
		cand_curve.geometry = geometry
		cand_curve.plane_point = plane_point
		cand_curve.plane_normal = plane_normal
		candidate_curves.append(cand_curve)
		creation_intersections.append([comb[0], comb[1], comb[2]])

	#print(len(candidate_curves))
	cand_curve_bbox = [tools.bbox_3d_single_stroke(cand_curve.geometry)
					   for cand_curve in candidate_curves]
	#print("collect_planes time: " + str(
	#	(time.clock() - start_time) / 60.0) + " min")
	#for cand_curve in candidate_curves:
	#	if np.all(np.isclose(cand_curve.plane_normal, 0.0)):
	#		print("cand_curve.plane_normal")
	#		print(cand_curve.plane_normal)
	# get intersection sets for all candidate lines: all intersections which are
	# within 0.1*length(cand_line)
	#start_time = time.clock()
	intersection_sets = []
	empty_intersection_sets = []
	for cand_curve_id, cand_curve in enumerate(candidate_curves):
		intersection_set = creation_intersections[cand_curve_id]
		line_length = tools_3d.line_3d_length(cand_curve.geometry)
		merge_dist = min(0.02 * bbox_diag, 0.1 * line_length)
		#merge_dist = min(0.005 * bbox_diag, 0.05 * line_length)
		for inter_3d_id, inter_3d in enumerate(intersections_3d):
			if inter_3d_id in intersection_set or inter_3d.coords_3d is None:
				continue
			if tools_3d.distance_to_bbox(inter_3d.coords_3d, cand_curve_bbox[cand_curve_id]) > merge_dist:
				continue
			#dist_old = tools_3d.distance_point_to_polyline(inter_3d.coords_3d,
			#										   cand_curve.geometry)

			dist = tools_3d.distance_point_to_polyline_vectorized(inter_3d.coords_3d,
																	   cand_curve.geometry)
			if np.isclose(dist, -1.0):
				continue
			#dist = tools_3d.distance_point_to_polyline(inter_3d.coords_3d,
			#										   cand_curve.geometry)
			#sys.exit()

			if dist < merge_dist:
			#if dist < 0.1 * line_length:
			#if dist < 0.02 * bbox_diag and dist < 0.1 * line_length:
				intersection_set.append(inter_3d_id)
		if len(intersection_set) > 0:
			intersection_sets.append(intersection_set)
		else:
			empty_intersection_sets.append(cand_curve_id)
	# remove empty candidate lines
	for del_id in sorted(empty_intersection_sets, reverse=True):
		del candidate_curves[del_id]

	#print("collect_intersections time: " + str(
	#	(time.clock() - start_time) / 60.0) + " min")

	#start_time = time.clock()
	clustered_candidate_curves, clustered_intersection_sets = \
		cluster_candidate_curves_v2(candidate_curves, intersection_sets, dep_node_id=s_id,
								 intersections_3d=intersections_3d)
	#print("clustering time: " + str(
	#	(time.clock() - start_time) / 60.0) + " min")

	#for cand_curve in clustered_candidate_curves:
	#	if np.all(np.isclose(cand_curve.plane_normal, 0.0)):
	#		print("clustered_cand_curve.plane_normal")
	#		print(cand_curve.plane_normal)
	return clustered_candidate_curves, clustered_intersection_sets

def cluster_candidate_curves_v2(candidate_curves, intersection_sets, dep_node_id=-1,
								intersections_3d=[]):
	#start_time = time.clock()
	clustered_candidate_curves = []
	clustered_intersection_sets = []
	potential_clusters = []
	intersection_sets_lengths = np.array([len(inter_set) for inter_set in intersection_sets])
	sorted_lenghts_ids = np.argsort(intersection_sets_lengths)
	_, meta_cluster_ids = np.unique(intersection_sets_lengths[sorted_lenghts_ids], return_index=True)
	meta_cluster_ids = np.append(meta_cluster_ids, len(intersection_sets_lengths))

	for vec_id in range(len(meta_cluster_ids)-1):
		#print(meta_cluster_ids[vec_id])
		#print(meta_cluster_ids[vec_id+1])
		#print(sorted_lenghts_ids[meta_cluster_ids[vec_id]:meta_cluster_ids[vec_id+1]])
		meta_cluster = [intersection_sets[i] for i in sorted_lenghts_ids
			[meta_cluster_ids[vec_id]:meta_cluster_ids[vec_id+1]].tolist()]
		meta_cluster = np.array(meta_cluster)
		meta_cluster = np.sort(meta_cluster, axis=1)
		unique_elems = np.unique(meta_cluster, axis=0)
		for unique_elem in unique_elems:
			sub_cluster_ids = np.where((meta_cluster == unique_elem).all(axis=1))
			potential_clusters.append(sorted_lenghts_ids
									  [meta_cluster_ids[vec_id]:meta_cluster_ids[vec_id+1]]
									  [sub_cluster_ids])

	checksum = 0
	all_ids = []
	for potential_cluster in potential_clusters:
		all_ids += list(potential_cluster)
		#print("potential_cluster: ", potential_cluster)
		#potential_cluster = [intersection_sets[i] for i in potential_cluster]
		#checksum += np.sum(potential_cluster)
		plane_normals = np.array([candidate_curves[i].plane_normal
								  for i in potential_cluster])
		#print(plane_normals)
		#print(plane_normals+1)
		#plane_normals = plane_normals + 1
		normal_adj_mat = (plane_normals[:, np.newaxis]*plane_normals).sum(axis=-1)
		normal_adj_mat = np.greater(np.deg2rad(5)/np.pi, 1.0 - np.abs(normal_adj_mat))
		#print(normal_adj_mat)
		#print(normal_adj_mat.shape)

		for cc in nx.connected_components(nx.Graph(normal_adj_mat)):
			cc_ids = np.array(potential_cluster)[list(cc)]
			#print(cc_ids)
			clustered_candidate_curves.append([candidate_curves[cc_id] for cc_id in cc_ids])
			clustered_intersection_sets.append(intersection_sets[potential_cluster[0]])
	#print("cluster merging time: " + str(
	#	(time.clock() - start_time) / 60.0) + " min")

	#print("end_result")
	#print(len(clustered_candidate_curves))
	#print(len(clustered_intersection_sets))

#	start_time = time.clock()
#	clustered_candidate_curves = []
#	clustered_intersection_sets = []
#	for inter_set_id, inter_set in enumerate(intersection_sets):
#		found_cluster = False
#		for inter_cluster_id, inter_cluster in enumerate(
#				clustered_intersection_sets):
#			#continue
#			if len(inter_set) != len(inter_cluster):
#				continue
#			if np.sum(np.in1d(inter_set, inter_cluster)) == len(inter_set):
#				plane_normal_1 = clustered_candidate_curves[inter_cluster_id][0].plane_normal
#				plane_normal_1 /= np.linalg.norm(plane_normal_1)
#				plane_normal_2 = candidate_curves[inter_set_id].plane_normal
#				plane_normal_2 /= np.linalg.norm(plane_normal_2)
#				if 1.0 - np.abs(np.dot(plane_normal_1, plane_normal_2)) < np.deg2rad(5)/np.pi:
#					found_cluster = True
#					clustered_candidate_curves[inter_cluster_id].append(
#						candidate_curves[inter_set_id])
#					break
#		if found_cluster:
#			continue
#		clustered_candidate_curves.append([candidate_curves[inter_set_id]])
#		clustered_intersection_sets.append(inter_set)
#
#	print("cluster intersection sets time: " + str(
#		(time.clock() - start_time) / 60.0) + " min")
#	start_time = time.clock()
	for cluster_id, clustered_curves in enumerate(clustered_candidate_curves):
		#if dep_node_id == 111:
		#	clustered_candidate_curves[cluster_id] = tools_3d.merge_n_curves(
		#		clustered_curves, VERBOSE=True, intersections=[intersections_3d[inter_id] for inter_id in clustered_intersection_sets[cluster_id]])
		#else:
		clustered_candidate_curves[cluster_id] = tools_3d.merge_n_curves(clustered_curves)
#	print("cluster merging time: " + str(
#		(time.clock() - start_time) / 60.0) + " min")
#
#	print("end_result")
#	print(len(clustered_candidate_curves))
#	print(len(clustered_intersection_sets))
#	sys.exit()
	return clustered_candidate_curves, clustered_intersection_sets

def cluster_candidate_curves(candidate_curves, intersection_sets, dep_node_id=-1,
							 intersections_3d=[]):
	start_time = time.clock()
	clustered_candidate_curves = []
	clustered_intersection_sets = []
	for inter_set_id, inter_set in enumerate(intersection_sets):
		found_cluster = False
		for inter_cluster_id, inter_cluster in enumerate(
				clustered_intersection_sets):
			#continue
			if len(inter_set) != len(inter_cluster):
				continue
			if np.sum(np.in1d(inter_set, inter_cluster)) == len(inter_set):
				plane_normal_1 = clustered_candidate_curves[inter_cluster_id][0].plane_normal
				plane_normal_1 /= np.linalg.norm(plane_normal_1)
				plane_normal_2 = candidate_curves[inter_set_id].plane_normal
				plane_normal_2 /= np.linalg.norm(plane_normal_2)
				if 1.0 - np.abs(np.dot(plane_normal_1, plane_normal_2)) < np.deg2rad(5)/np.pi:
					found_cluster = True
					clustered_candidate_curves[inter_cluster_id].append(
						candidate_curves[inter_set_id])
					break
		if found_cluster:
			continue
		clustered_candidate_curves.append([candidate_curves[inter_set_id]])
		clustered_intersection_sets.append(inter_set)

	print("cluster intersection sets time: " + str(
		(time.clock() - start_time) / 60.0) + " min")
	start_time = time.clock()
	for cluster_id, clustered_curves in enumerate(clustered_candidate_curves):
		#if dep_node_id == 111:
		#	clustered_candidate_curves[cluster_id] = tools_3d.merge_n_curves(
		#		clustered_curves, VERBOSE=True, intersections=[intersections_3d[inter_id] for inter_id in clustered_intersection_sets[cluster_id]])
		#else:
		clustered_candidate_curves[cluster_id] = tools_3d.merge_n_curves(clustered_curves)
	print("cluster merging time: " + str(
		(time.clock() - start_time) / 60.0) + " min")

	return clustered_candidate_curves, clustered_intersection_sets

# all_context_lines: [ [geom, dep_inter_set, sketch_versions],
# 						...
# 						...									]
# merged_context_lines: [ [geom, [dep_inter_set_1, ...], [sketch_versions_1, ...]],
#						...
#						...									]
# merge geometries if endpoints within 10% of length of shorter stroke
def merge_over_contexts_curves(all_context_curves):

	# return identity for now
	merged_context_curves = []
	for context_curve_id, context_curve in enumerate(all_context_curves):
		merged_context_curves.append([context_curve[0], [context_curve[1]],
									 [context_curve[2]]])
	return merged_context_curves

def merge_over_contexts_straight_lines(all_context_lines, bbox_diag, dep_node_id=-1):
	for c in all_context_lines:
		if len(c[-1]) == 0:
			print("beginning: EMPTY_VERSION SET")
	merged_context_lines = []
	clustered_geom_lines = []
	for context_line_id, context_line in enumerate(all_context_lines):
		found_cluster = False
		for cluster_id, cluster in enumerate(clustered_geom_lines):
			for geom_line_id in cluster:
				shorter_length = min(
					tools_3d.line_3d_length(all_context_lines[geom_line_id][0].geometry),
					tools_3d.line_3d_length(
						all_context_lines[context_line_id][0].geometry))
				merge_dist = min(0.005*bbox_diag, 0.05*shorter_length)
				merge_dist = min(0.02*bbox_diag, 0.1*shorter_length)
				if (np.linalg.norm(all_context_lines[geom_line_id][0].geometry[0] -
									all_context_lines[context_line_id][0].geometry[
										0]) < merge_dist) and \
					(np.linalg.norm(all_context_lines[geom_line_id][0].geometry[-1] -
									all_context_lines[context_line_id][0].geometry[
										-1]) < merge_dist):#) or \
						#((np.linalg.norm(all_context_lines[geom_line_id][0].geometry[0] -
						#				 all_context_lines[context_line_id][0].geometry[
						#					 -1]) < 0.1 * shorter_length) and
						# (np.linalg.norm(
						#	 all_context_lines[geom_line_id][0].geometry[-1] -
						#	 all_context_lines[context_line_id][0].geometry[
						#		 0]) < 0.1 * shorter_length)):
					plane_normal_1 = all_context_lines[geom_line_id][0].plane_normal
					plane_normal_2 = all_context_lines[context_line_id][0].plane_normal
					if np.isclose(np.linalg.norm(plane_normal_1), 0.0) or \
							np.isclose(np.linalg.norm(plane_normal_2), 0.0):
						continue
					plane_normal_1 /= np.linalg.norm(plane_normal_1)
					plane_normal_2 /= np.linalg.norm(plane_normal_2)
					if 1.0 - np.abs(np.dot(plane_normal_1,
										   plane_normal_2)) < np.deg2rad(5) / np.pi:
						found_cluster = True
						clustered_geom_lines[cluster_id].append(context_line_id)
						break
			if found_cluster:
				break
		if not found_cluster:
			clustered_geom_lines.append([context_line_id])

	#print("len(clustered_geom_lines)")
	#print(len(clustered_geom_lines))
	for cluster_id, cluster in enumerate(clustered_geom_lines):
		#print("len(cluster): ", len(cluster))
		new_cluster = []
		clustered_curves = [all_context_lines[geom_line_id][0] for geom_line_id in cluster]
		#clustered_intersections = []
		#for inter_set in [all_context_lines[geom_line_id][1] for geom_line_id in cluster]:
		#	for inter in inter_set:
		#		clustered_intersections.append(inter)
		#geom = tools_3d.merge_n_curves(
		#	[all_context_lines[geom_line_id][0] for geom_line_id in cluster])
		#if dep_node_id == 42:
		#	geom = tools_3d.merge_n_curves(
		#			clustered_curves, VERBOSE=True, intersections=clustered_intersections)
		#else:
		geom = tools_3d.merge_n_curves(clustered_curves)
		new_cluster.append(geom)
		new_cluster.append(
			[all_context_lines[geom_line_id][1] for geom_line_id in cluster])
		new_cluster.append(
			[all_context_lines[geom_line_id][2] for geom_line_id in cluster])
		merged_context_lines.append(new_cluster)

	# cluster identical sketch_version_sets
	for context_lined_id, (_, inter_sets, version_sets) in enumerate(
			merged_context_lines):
		version_set_clusters = []
		inter_set_clusters = []
		cluster_ids = []
		for version_set_id, version_set in enumerate(version_sets):
			if len(version_set) == 0:
				print("EMPTY VERSION SET!")
			found_cluster = False
			for version_set_cluster_id, version_set_cluster in enumerate(
					version_set_clusters):
				if len(version_set_cluster) == len(version_set) and \
						np.sum(
							np.in1d(version_set_cluster, version_set)) == len(
					version_set):
					found_cluster = True
					cluster_ids[version_set_cluster_id].append(version_set_id)
					inter_sets[version_set_cluster_id] += inter_sets[
						version_set_id]
				if found_cluster:
					break
			if not found_cluster:
				version_set_clusters.append(version_set)
				inter_set_clusters.append(inter_sets[version_set_id])
				cluster_ids.append([version_set_id])

			# kick out identical intersections
			for inter_set_cluster_id, inter_set_cluster in enumerate(
					inter_set_clusters):
				inter_ids = [inter.inter_id for inter in inter_set_cluster]
				u, inter_ids_unique = np.unique(inter_ids, return_index=True)
				inter_set_clusters[inter_set_cluster_id] = [
					inter_set_cluster[inter_id] for inter_id in
					inter_ids_unique]
			merged_context_lines[context_lined_id][1] = inter_set_clusters
			merged_context_lines[context_lined_id][2] = version_set_clusters
			if len(version_set_clusters) == 0:
				print("EMPTY VERSION CLUSTER!")

	# merged_context_lines = []
	return merged_context_lines

def project_tangents_to_3d(context_curves, dep_graph, s_id, camera):
	for context_curve in context_curves:
		curve = context_curve[0]
		for inter_id, inter in enumerate(context_curve[1]):
			curr_s_id = np.argwhere(np.array(inter.stroke_ids) == s_id).flatten()[0]
			curve_tangent = inter.tangents_2d[1]
			curve_tangent_3d = camera.lift_polyline_to_plane_vectorized(curve_tangent,
															 curve.plane_point,
															 curve.plane_normal)

			other_s_id = inter.stroke_ids[1-curr_s_id]
			other_tangent = inter.tangents_2d[0]
			sketch_version = context_curve[2][0]
			other_dep_node_id = dep_graph.stroke_id_to_dep_node_id[other_s_id]
			other_cand_node_id = dep_graph.sketch_versions_reference \
				[sketch_version][other_dep_node_id][0]
			other_line = dep_graph.dependency_nodes[other_dep_node_id].candidate_nodes[other_cand_node_id]
			if dep_graph.dependency_nodes[dep_graph.stroke_id_to_dep_node_id[other_s_id]].is_curve:
				other_tangent_3d = camera.lift_polyline_to_plane_vectorized(other_tangent,
																 other_line.plane_point,
																 other_line.plane_normal)
			else:
				other_tangent_3d = other_line.geometry

			curve_tangent_3d = np.array(curve_tangent_3d[-1] - curve_tangent_3d[0])
			other_tangent_3d = np.array(other_tangent_3d[-1] - other_tangent_3d[0])
			curve_tangent_3d /= np.linalg.norm(curve_tangent_3d)
			other_tangent_3d /= np.linalg.norm(other_tangent_3d)
			angle = 180*acos(min(1.0, abs(np.dot(other_tangent_3d, curve_tangent_3d))))/np.pi
			context_curve[1][inter_id].tangents_3d = [other_tangent_3d, curve_tangent_3d]
			context_curve[1][inter_id].tangents_angle_3d = angle


def filter_out_foreshortened_lines(candidate_curves, intersections_sets,
								   dep_graph, s_id, strokes_topology, median_length_2d,
								   bbox):
	delete_indices = []
	local_candidate_lines = deepcopy(candidate_curves)
	local_intersections_sets = deepcopy(intersections_sets)
	for cand_line_id, cand_line in enumerate(candidate_curves):
		#print(cand_line.geometry < bbox[:3])
		#print(np.any(cand_line.geometry < bbox[:3]))
		#print(np.any(cand_line.geometry < bbox[:3]) or np.any(cand_line.geometry > bbox[3:]))
		if tools_3d.line_3d_length(cand_line.geometry) / dep_graph.median_length_3d > \
				2.0 * strokes_topology[s_id]["linestring"].length / median_length_2d or \
				np.any(cand_line.geometry < bbox[:3]) or np.any(cand_line.geometry > bbox[3:]):
			delete_indices.append(cand_line_id)
		#if np.any(cand_line.geometry < bbox[:3]) or np.any(cand_line.geometry > bbox[3:]):
		#	delete_indices.append(cand_line_id)
	if len(delete_indices) == len(candidate_curves):
		return local_candidate_lines, local_intersections_sets
	for del_id in sorted(delete_indices, reverse=True):
		del local_candidate_lines[del_id]
		del local_intersections_sets[del_id]
	return local_candidate_lines, local_intersections_sets
