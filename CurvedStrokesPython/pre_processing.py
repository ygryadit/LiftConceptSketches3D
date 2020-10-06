import tools
import os
import numpy as np
import compute_intersections
import reproject_2d_on_3d_linestring
from shapely.geometry import LineString

ACC_RADIUS_SCALE_FACTOR = 2.0
TANGENT_SENSIBILITY = 0.4
BBOX_SCALE_FACTOR = 0.10

def pre_process(file_name, designer_name, object_name, strokes_topology,
				intersections, cam_params):
	# PRE-PROCESSING
	tools.setup_dict_keys_strokes_topology(strokes_topology)
	strokes_array = tools.convert_strokes_topology_to_strokes_array(
		strokes_topology)
	strokes_array_3d = tools.convert_strokes_topology_to_strokes_array(strokes_topology, 1)
	lines_array = tools.convert_strokes_array_to_lines_array(strokes_array)
	lines_array_3d = tools.convert_strokes_array_to_lines_array(strokes_array_3d)
	# fill in strokes_array_3d and lines_array_3d in strokes_topology
	for s_idx, s in enumerate(strokes_topology):
		strokes_topology[s_idx]["stroke"] = strokes_array[s_idx]
		strokes_topology[s_idx]["linestring"] = lines_array[s_idx]
		strokes_topology[s_idx]["stroke_3d"] = strokes_array_3d[s_idx]
		strokes_topology[s_idx]["linestring_3d"] = lines_array_3d[s_idx]
		strokes_topology[s_idx]["curve_intersections"] = []

	cluster_folder = os.path.join("full_clusters/", designer_name, object_name,
								  "view1", "lines_separation")
	cluster_file_name = os.path.join(cluster_folder, "curves.cluster")
	if os.path.exists(cluster_file_name):
		cluster_geometry_file_name = os.path.join(cluster_folder,
												  "curves_fit.scap")
		offset_file_name = os.path.join(cluster_folder, "curves.offset")
	else:
		print(cluster_file_name + " does not exist! ")
		USE_CLUSTERS = False
	strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology)
	tools.fit_ellipses(strokes_topology, strokes_array)
	strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology)
	if cluster_file_name != "toto.cluster":
		tools.scale_acc_radius(strokes_topology, ACC_RADIUS_SCALE_FACTOR)
		clusters, cluster_scap_ids, cluster_lines = tools.read_stroke_clusters(
			cluster_file_name, cluster_geometry_file_name, offset_file_name)

		del_clusters = []
		for cluster_idx, cluster in enumerate(clusters):
			if len(cluster) < 2:
				# we don't want to use aggregated curves when they are just
				# replacing a single input curves
				del_clusters.append(cluster_idx)
			for c in cluster:
				if strokes_topology[c]["is_ellipse"]:
					del_clusters.append(cluster_idx)
					break

		non_planar_clusters = []
		# check which clusters are non-planar and break them apart
		check_non_planar_indices = []
		for cluster_idx, cluster in enumerate(clusters):
			# get all tangential intersection with the scaffold lines
			if len(cluster) < 2 or cluster_idx in del_clusters:
				continue
			cluster = np.sort(cluster)
			curve_idx = cluster[-1]
			acc_radius = strokes_topology[curve_idx]["accuracy_radius"]
			_, is_ellipse = tools.check_is_ellipse(np.array(cluster_lines[cluster_idx]),
												   acc_radius)
			if is_ellipse:
				# stroke-aggregator doesn't aggregate ellipses that well
				# todo: actually, we should check if any curve in the cluster
				# can be identified as an ellipse
				del_clusters.append(cluster_idx)
				continue

			check_non_planar_indices.append(curve_idx)
			curve_intersections = compute_intersections.compute_curve_2D_intersections_single_curve_input_linestring_simplified(
				curve_idx, LineString(cluster_lines[cluster_idx]), acc_radius,
				strokes_array, lines_array, lines_array_3d,
				strokes_topology, tangent_sensibility=TANGENT_SENSIBILITY)
			# get 3D tangents
			curve_intersections = reproject_2d_on_3d_linestring. \
				get_3d_constraints_single_curve(curve_intersections, [],
												strokes_array_3d, [],
												cam_params, cc_intersections=False)
			for c in cluster:
				strokes_topology[c]["curve_intersections"] = curve_intersections
				strokes_topology[c]["stroke"] = cluster_lines[cluster_idx]
				strokes_topology[c]["linestring"] = LineString(cluster_lines[cluster_idx])
		tools.detect_non_planar_candidates(strokes_topology,
										   curve_indices=check_non_planar_indices)

		# check which clusters are non-planar
		for cluster_idx, cluster in enumerate(clusters):
			if len(cluster) > 1 and strokes_topology[cluster[-1]]["non_planar_curve"]:
				del_clusters.append(cluster_idx)
				non_planar_clusters.append(cluster_idx)
		del_clusters = np.unique(del_clusters).tolist()

		for del_idx in sorted(del_clusters, reverse=True):
			del clusters[del_idx]
			del cluster_scap_ids[del_idx]
			del cluster_lines[del_idx]

		strokes_array = tools.convert_strokes_topology_to_strokes_array(
			strokes_topology)
		old_new_stroke_mapping = tools.regroup_strokes(strokes_topology,
							  	clusters, cluster_scap_ids, cluster_lines,
							  	intersections)
		tools.scale_acc_radius(strokes_topology, 1.0/ACC_RADIUS_SCALE_FACTOR)

		strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology)

	tools.clean_strokes_topology(strokes_topology, intersections, only_markers=True)
	strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology)

	tools.remove_markers(strokes_topology, intersections)
	strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology)


	cam_params["view_dir"] = np.array(cam_params["view_dir"])
	strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology)
	tools.fill_in_cam_params_inv_mats(cam_params)
	cam_params["im_dim"] = tools.get_im_dim(strokes_array)
	tools.recompute_im_vec(cam_params)

	# 2D cubic bezier representations for curve intersections
	tools.fit_beziers(strokes_topology, strokes_array, VERBOSE=0)

	# clustering of easily identifiable oversketched curves
	lines_array = tools.convert_strokes_array_to_lines_array(strokes_array)
	strokes_array_3d = tools.convert_strokes_topology_to_strokes_array(
		strokes_topology, 1)
	lines_array_3d = tools.convert_strokes_array_to_lines_array(
		strokes_array_3d)
	tools.scale_acc_radius(strokes_topology, ACC_RADIUS_SCALE_FACTOR)
	# fill strokes_topology with information
	for s_idx, s in enumerate(strokes_topology):
		if s["primitive_type"] != 1:
			continue

		curve_intersections = compute_intersections.compute_curve_2D_intersections_single_curve(
			s_idx, strokes_array, strokes_array_3d,
			lines_array, lines_array_3d, strokes_topology, intersections,
			cc_intersections=True, tangent_sensibility=TANGENT_SENSIBILITY,
			reject=False)
		strokes_topology[s_idx]["curve_intersections"] = curve_intersections

	parent_list, child_list = tools.get_absorbed_strokes(strokes_topology)
	for child_idx in sorted(child_list, reverse=True):
		parent_id = np.argwhere(np.array(child_list) == child_idx).flatten()[0]
		old_new_stroke_mapping.append([strokes_topology[child_idx]["old_index"],
									   strokes_topology[parent_id]["old_index"]])
		del strokes_topology[child_idx]
	strokes_array = tools.convert_strokes_topology_to_strokes_array(
		strokes_topology)
	tools.scale_acc_radius(strokes_topology, 1.0 / ACC_RADIUS_SCALE_FACTOR)

	tools.fit_ellipses(strokes_topology, strokes_array)
	strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology)
	tools.scale_acc_radius(strokes_topology, ACC_RADIUS_SCALE_FACTOR)

	# more convenient to have directly acces to the 2D coordinates
	strokes_array_3d = tools.convert_strokes_topology_to_strokes_array(strokes_topology, 1)
	lines_array = tools.convert_strokes_array_to_lines_array(strokes_array)
	lines_array_3d = tools.convert_strokes_array_to_lines_array(strokes_array_3d)
	tools.update_strokes_indices(intersections, strokes_topology)
	tools.get_stroke_planes(strokes_topology, intersections, strokes_array_3d)
	bbox = tools.intersections_bbox_3d(intersections)
	tools.scale_bbox(bbox, BBOX_SCALE_FACTOR)
	strokes_array = tools.convert_strokes_topology_to_strokes_array(strokes_topology)
	# fill in strokes_array_3d and lines_array_3d in strokes_topology
	for s_idx, s in enumerate(strokes_topology):
		strokes_topology[s_idx]["stroke"] = strokes_array[s_idx]
		strokes_topology[s_idx]["linestring"] = lines_array[s_idx]
		strokes_topology[s_idx]["stroke_3d"] = strokes_array_3d[s_idx]
		strokes_topology[s_idx]["linestring_3d"] = lines_array_3d[s_idx]

	# pre-compute curve_intersections ONLY with the scaffold and put them into
	# strokes_topology
	# get curve indices
	curve_indices = []
	for stroke_idx, stroke in enumerate(strokes_topology):
		if stroke["primitive_type"] == 1:
			curve_indices.append(stroke_idx)
	curve_indices = np.array(curve_indices)
	CC_INTERSECTIONS = True
	precompute_curves_3d = []
	del_curve_indices = []
	empty_curve_indices = []
	for c_idx_idx, c_idx in enumerate(curve_indices):

		strokes_topology[c_idx]["curve_intersections"] = []
		# get 2D intersections
		curve_intersections = compute_intersections.compute_curve_2D_intersections_single_curve(
			c_idx, strokes_array, strokes_array_3d,
			lines_array, lines_array_3d, strokes_topology, intersections,
			cc_intersections=True, tangent_sensibility=TANGENT_SENSIBILITY)
		# filter out intersections with empty curves
		for i in reversed(range(len(curve_intersections))):
			if curve_intersections[i]["straight_line_idx"] in empty_curve_indices:
				del curve_intersections[i]
		# also filter out intersections that are too much overlapping with other
		# curves. This is to avoid a "being glued together" effect
		tools.filter_out_curve_intersections(curve_intersections,
											 strokes_topology)

		if len(curve_intersections) == 0:
			del_curve_indices.append(c_idx_idx)
			empty_curve_indices.append(c_idx)
			continue
		# get 3D intersections for scaffold intersections only
		curve_intersections = reproject_2d_on_3d_linestring. \
			get_3d_constraints_single_curve(curve_intersections, strokes_array,
											strokes_array_3d, strokes_topology,
											cam_params, cc_intersections=False,
											only_cc_intersections=False)
		if len(curve_intersections) == 0:
			del_curve_indices.append(c_idx_idx)
			empty_curve_indices.append(c_idx)
			continue
		strokes_topology[c_idx]["curve_intersections"] = curve_intersections

	curve_indices = np.delete(curve_indices, del_curve_indices)

	# get curve cluster based on simple intersections
	# those clusters will each be treated by a separate tree
	clusters = tools.cluster_intersecting_strokes_networkx_v2(
		curve_indices, strokes_topology)

	stroke_transformations = []
	for old_new_s in old_new_stroke_mapping:
		if len(stroke_transformations) > 0:
			already_existing_indices = np.argwhere(np.array(stroke_transformations)[:, 1] == old_new_s[0]).flatten()
			if len(already_existing_indices) > 0:
				for alread_i in already_existing_indices:
					stroke_transformations[alread_i][1] = old_new_s[0]
				continue
		stroke_transformations.append(old_new_s)

	return clusters, bbox, np.array(stroke_transformations)
