import numpy as np
import shapely
from shapely.geometry import LineString
from shapely.ops import nearest_points
import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+"/fitCurves/")
import fitCurves
import bezier as bz
from scipy.special import binom
from rdp import rdp
from interval import interval
import networkx as nx

# fetches only the 2D information
def convert_strokes_topology_to_strokes_array(strokes_topology, three_d=0, use_rdp=True):
    strokes_array = []
    for stroke_idx, stroke in enumerate(strokes_topology):
        stroke_array = []
        if three_d:
            if not stroke["depth_assigned"]:
                strokes_array.append(np.array(stroke_array))
                continue
            if use_rdp and stroke["primitive_type"] == 0:
                p0 = np.array(stroke["primitive_geom_3D"][0])
                p1 = np.array(stroke["primitive_geom_3D"][1])
                stroke_array = [p0, p1]
            else:
                for p_idx, p in enumerate(stroke["points3D"]):
                    stroke_array.append(np.array([p[0], p[1], p[2]]))
        else:
            if stroke["primitive_type"] == 0:
                p0 = np.array([stroke["primitive_geom"][0], stroke["primitive_geom"][2]])
                p1 = np.array([stroke["primitive_geom"][1], stroke["primitive_geom"][3]])
                stroke_array = [p0, p1]
            else:
                if isinstance(stroke["points2D"], dict):
                    stroke_array.append([])
                elif use_rdp and stroke["is_ellipse"]:
                    stroke_array = np.array(stroke["bezier_linestring"])
                else:
                    for p in stroke["points2D"]:
                        if isinstance(p, dict):
                            stroke_array.append(np.array([p["x"], p["y"]]))
                        else:
                            stroke_array.append(np.array(p))
                if use_rdp:
                    stroke_array = rdp(stroke_array, epsilon=1.0)
        strokes_array.append(np.array(stroke_array))

    return strokes_array

def convert_strokes_array_to_lines_array(strokes_array):
    lines_array = []
    for stroke in strokes_array:
        if len(stroke) < 2:
            lines_array.append(LineString())
        else:
            lines_array.append(LineString(stroke))
    return lines_array

def get_intersection_list(line_a, line_b, radius):
    multi_point = line_a.buffer(radius).intersection(line_b)
    intersection_point_list = []
    # the intersection can have different types, but we are only
    # interested in the coordinates
    if isinstance(multi_point, shapely.geometry.MultiLineString):
        for line_string in list(multi_point.geoms):
            intersection_point_list.append(np.array(line_string))
    else:
        intersection_point_list.append(np.array(multi_point))
    return intersection_point_list

def is_tangential(l1, l2, acceptance_rate=0.3, acc_radius=-1.0, angle=0.95):

    step_size = 0.05
    tangential_score = 0
    nb_probes = 0
    # go through l1 and compare tangents with nearest point on l2
    t = 0.0

    while t < 1.0:

        p = l1.interpolate(t, normalized=True)
        p_prim = l1.interpolate(t+step_size, normalized=True)
        p_coords = np.array(p)
        p_prim_coords = np.array(p_prim)

        t += step_size
        p, q = nearest_points(p, l2)
        p_prim, q_prim = nearest_points(p_prim, l2)

        q_coords = np.array(q)
        q_prim_coords = np.array(q_prim)

        if np.all(np.isclose(p,p_prim)) or np.all(np.isclose(q,q_prim)):
            continue
        t1 = p_prim_coords - p_coords
        t2 = q_prim_coords - q_coords
        cos_alpha = np.dot(t1, t2)/(np.linalg.norm(t1)*np.linalg.norm(t2))
        if cos_alpha > angle:
            tangential_score += 1
        nb_probes += 1

    if nb_probes > 1:
        return float(tangential_score)/float(nb_probes) >= acceptance_rate
    else:
        return False

# tests for tangentiality on extended versions of l2
# this is basically a test for "tangential continuity"
def is_tangential_extended(l1, l2, acceptance_rate=0.3):
    # check which end of l2 is nearest to l1
    l2_start_point = l2.interpolate(0.0, normalized=True)
    l2_end_point = l2.interpolate(1.0, normalized=True)
    # double l2's length in nearest point direction
    l2_closest_point = l2_start_point
    l2_farthest_point = l2_end_point
    if l1.distance(l2_end_point) < l1.distance(l2_start_point):
        l2_closest_point = l2_end_point
        l2_farthest_point = l2_start_point
    extension_vec = np.array(l2_closest_point) - np.array(l2_farthest_point)
    l2_extended = shapely.geometry.LineString([extension_vec + l2_closest_point,
                                              l2_farthest_point])
    return is_tangential(l1, l2_extended, acceptance_rate)

def update_strokes_indices(intersections, strokes_topology):
    for stroke_idx, stroke in enumerate(strokes_topology):
        for inter_idx, inter in enumerate(intersections):
            if stroke["old_index"] == inter["strokes_indices"][0] - 1:
                intersections[inter_idx]["strokes_indices"][0] = stroke_idx
            if stroke["old_index"] == inter["strokes_indices"][1] - 1:
                intersections[inter_idx]["strokes_indices"][1] = stroke_idx

# intersected_line is a LineString
def is_trihedral_end(stroke_idx, intersected_line_part_curved,
                     intersected_line_part_straight, intersections,
                     acc_radius, strokes_array, strokes_topology):
    # check if at least one other non-colinear intersection is close-by
    valence_count = 0
    for straight_inter_idx, straight_inter in enumerate(intersections):
        straight_inter_p = shapely.geometry.Point(straight_inter["coordinates2D"])
        if intersected_line_part_curved.distance(straight_inter_p) < acc_radius\
                and straight_inter["is_active"]:
            if not (straight_inter["strokes_indices"][0] < stroke_idx and \
                    straight_inter["strokes_indices"][1] < stroke_idx):
                continue
#            # to get the old straight lines, we have to look for the old indices
            index_1, index_2 = straight_inter["strokes_indices"]
            l1 = strokes_array[index_1]
            l2 = strokes_array[index_2]
            l1_linestring = shapely.geometry.LineString(l1)
            l2_linestring = shapely.geometry.LineString(l2)
            if not(l1_linestring.buffer(1.0).contains(intersected_line_part_straight.buffer(0.9)) or \
                   l2_linestring.buffer(1.0).contains(intersected_line_part_straight.buffer(0.9))):
                continue
            if not straight_inter["collinear"]:
                valence_count += 1
    if valence_count > 0:
        return True

    return False

def stroke_start_covered(stroke_idx, stroke_part, strokes_array, acc_radius):
    stroke_linestring = shapely.geometry.LineString(strokes_array[stroke_idx])
    stroke_start = substring(stroke_linestring, 0.0, 0.10, normalized=True)

    return stroke_start.distance(stroke_part) < 1.0 * acc_radius

def stroke_end_covered(stroke_idx, stroke_part, strokes_array, acc_radius):
    stroke_linestring = shapely.geometry.LineString(strokes_array[stroke_idx])
    stroke_end = substring(stroke_linestring, 0.90, 1.0, normalized=True)

    return stroke_end.distance(stroke_part) < 1.0 * acc_radius

def project_3d_stroke(stroke, cam_params):
    proj = np.array(cam_params["P"])
    stroke_projected = []
    for p in stroke:
        v = np.append(np.array(p), 1.0)
        v_proj = np.dot(proj, v)
        v_proj_norm = v_proj[:2]/v_proj[-1]
        stroke_projected.append(v_proj_norm)
    stroke_projected = np.array(stroke_projected)
    return stroke_projected

#https://github.com/Toblerity/Shapely/blob/b2c562c6a78cf86fa284ca045946b389dfdf1fe8/shapely/ops.py
def substring(geom, start_dist_real, end_dist_real, normalized=False):
    """Return a line segment between specified distances along a linear geometry.
    Negative distance values are taken as measured in the reverse
    direction from the end of the geometry. Out-of-range index
    values are handled by clamping them to the valid range of values.
    If the start distances equals the end distance, a point is being returned.
    If the normalized arg is True, the distance will be interpreted as a
    fraction of the geometry's length.
    """

    start_dist = min(start_dist_real, end_dist_real)
    end_dist = max(start_dist_real, end_dist_real)
    assert (isinstance(geom, LineString))

    start_point = geom.interpolate(start_dist, normalized=normalized)
    end_point = geom.interpolate(end_dist, normalized=normalized)

    min_dist = min(start_dist, end_dist)
    max_dist = max(start_dist, end_dist)
    if normalized:
        min_dist *= geom.length
        max_dist *= geom.length

    vertex_list = [np.array(start_point)]
    coords = list(geom.coords)
    for p in coords:
        pd = geom.project(shapely.geometry.Point(p))
        if pd <= min_dist:
            pass
        elif min_dist < pd < max_dist:
            vertex_list.append(np.array(p))
        else:
            break
    vertex_list.append(np.array(end_point))

    # reverse direction of section
    if start_dist > end_dist:
        vertex_list = reversed(vertex_list)

    return LineString(vertex_list)

def get_3d_line_ends_in_correct_order(start_point, straight_3d_line, cam_params):
    start_point_3d = line_interpolate(straight_3d_line, 0.0)
    end_point_3d = line_interpolate(straight_3d_line, 1.0)
    start_end_projected = project_3d_stroke(np.array([np.array(start_point_3d),
                                                      np.array(end_point_3d)]),
                                            cam_params)
    if start_point.distance(shapely.geometry.Point(start_end_projected[0])) < \
            start_point.distance(
                shapely.geometry.Point(start_end_projected[1])):
        return start_point_3d, end_point_3d
    else:
        return end_point_3d, start_point_3d

def bezier_interpolation(endpoints, steps):
    interpolated_curve = []
    t = 0.0
    for i in range(steps+1):
        t = float(i)/float(steps)
        # cubic Bezier curve
        h1 = -1.0*t*t*t + 3.0*t*t - 3.0*t + 1
        h2 = 3.0*t*t*t - 6.0*t*t + 3.0*t
        h3 = -3.0*t*t*t + 3.0*t*t
        h4 = t*t*t
        p = h1*endpoints[0] + \
            h2 * endpoints[1] + \
            h3 * endpoints[2] + \
            h4 * endpoints[3]
        interpolated_curve.append(p)
        t += 0.1
    return interpolated_curve

def fit_beziers(strokes_topology, strokes_array, VERBOSE):

    for curve_idx, curve in enumerate(strokes_topology):
        if strokes_topology[curve_idx]["primitive_type"] != 1:
            continue
        if strokes_topology[curve_idx]["is_ellipse"] == 1:
            continue
        points = strokes_array[curve_idx]
        bezier_curves = fitCurves.fitCurve(points, strokes_topology[curve_idx]["accuracy_radius"])
        strokes_topology[curve_idx]["bezier_cps"] = bezier_curves


def fill_in_cam_params_inv_mats(cam_params, VERBOSE=0):

    if not "R_t_inv" in cam_params:
        r = np.array(cam_params["R"])
        t = np.array(cam_params["t"])
        r_t_inv = np.ones(shape=(4, 4))
        r_t_inv[:3, :3] = r
        r_t_inv[:3, 3] = t
        r_t_inv[3, :3] = 0.0
        r_t_inv = np.linalg.inv(r_t_inv)

        f = np.array(cam_params["f"])
        u0 = np.array(cam_params["principal_point"][0])
        v0 = np.array(cam_params["principal_point"][1])

        k = np.array([[f, 0.0, u0],
                      [0.0, f, v0],
                      [0.0, 0.0, 1.0]])
        k_inv = np.linalg.inv(k)
        cam_params["R_t_inv"] = r_t_inv
        cam_params["K_inv"] = k_inv

# expected value between [0.0, 1.0]
def line_interpolate(l, t, normalized=True):
    p0 = np.array(l)[0]
    p1 = np.array(l)[-1]
    p_interp = p0 + t*(p1-p0)
    return shapely.geometry.Point(p_interp)

def scale_bbox(bbox, scale_factor):

    bbox_min = bbox[:3]
    bbox_max = bbox[3:]
    bbox_mid = (bbox_min + bbox_max)/2.0
    bbox_min = bbox_min + (bbox_min-bbox_mid)*scale_factor
    bbox_max = bbox_max + (bbox_max-bbox_mid)*scale_factor
    bbox[:3] = bbox_min
    bbox[3:] = bbox_max

def scale_acc_radius(strokes_topology, acc_scale_factor):
    for s_idx, s in enumerate(strokes_topology):
        if s["accuracy_radius"] is not None:
            strokes_topology[s_idx]["accuracy_radius"] *= acc_scale_factor

def bbox_3d(strokes_array_3d):
    # min_x, min_y, min_z, max_x, max_y, max_z
    strokes_array_3d = np.array(strokes_array_3d)
    bbox = np.zeros(6, dtype=float)
    bbox[:3] = 1000.0
    bbox[3:] = -1000.0
    for s in strokes_array_3d:
        if len(s) == 0:
            continue
        bbox[0] = np.minimum(bbox[0], np.min(s[:, 0]))
        bbox[1] = np.minimum(bbox[1], np.min(s[:, 1]))
        bbox[2] = np.minimum(bbox[2], np.min(s[:, 2]))
        bbox[3] = np.maximum(bbox[3], np.max(s[:, 0]))
        bbox[4] = np.maximum(bbox[4], np.max(s[:, 1]))
        bbox[5] = np.maximum(bbox[5], np.max(s[:, 2]))
    return bbox

def bbox_3d_single_stroke(stroke_3d):
    # min_x, min_y, min_z, max_x, max_y, max_z
    stroke_3d = np.array(stroke_3d)
    bbox = np.zeros(6, dtype=float)
    bbox[:3] = 1000.0
    bbox[3:] = -1000.0
    bbox[0] = np.minimum(bbox[0], np.min(stroke_3d[:, 0]))
    bbox[1] = np.minimum(bbox[1], np.min(stroke_3d[:, 1]))
    bbox[2] = np.minimum(bbox[2], np.min(stroke_3d[:, 2]))
    bbox[3] = np.maximum(bbox[3], np.max(stroke_3d[:, 0]))
    bbox[4] = np.maximum(bbox[4], np.max(stroke_3d[:, 1]))
    bbox[5] = np.maximum(bbox[5], np.max(stroke_3d[:, 2]))
    return bbox

def intersections_bbox_3d(intersections):
    inter_array = []
    for inter in intersections:
        if inter["is_active"]:
            inter_array.append(np.array([inter["coordinates3D"]]))
    inter_array = np.array(inter_array)
    return bbox_3d(inter_array)

def get_closest_param_bezier(cps, p):
    stepsize = 0.01
    t_potentials = np.arange(0.0, 1.0+stepsize, stepsize)
    distances = np.zeros(len(t_potentials))
    for i in range(len(distances)):
        t_i = t_potentials[i]
        b_t_i = 0.0
        for k in range(4):
            b_k = binom(3, k) * (t_i ** k) * (1.0 - t_i) ** (3 - k)
            b_k *= cps[k]
            b_t_i += b_k
        distances[i] = np.linalg.norm(p - b_t_i)
    return t_potentials[np.argmin(distances)], np.min(distances)

def get_closest_param_beziers(beziers, p):

    if len(beziers) == 1:
        t, _ = get_closest_param_bezier(beziers[0], p)
        return t

    # get length of each curve
    lengths = []
    linestrings = []
    for bezier in beziers:
        linestring = LineString(bezier_interpolation(bezier, 100))
        linestrings.append(linestring)
        lengths.append(linestring.length)
    lengths = np.array(lengths)

    # compute total length
    total_length = np.sum(lengths)

    # get index of closest bezier part
    t_params = []
    distances = []
    for bezier in beziers:
        t, dist = get_closest_param_bezier(bezier, p)
        t_params.append(t)
        distances.append(dist)
    distances = np.array(distances)
    bezier_idx = np.argmin(distances)

    # compute final_t_param
    previous_length = np.sum(lengths[:bezier_idx])
    current_length = t_params[bezier_idx]*lengths[bezier_idx]
    final_t_param = (previous_length+current_length)/total_length
    return final_t_param

# remove all primitive_type -1
def setup_dict_keys_strokes_topology(strokes_topology):
    for i in range(len(strokes_topology)):
        strokes_topology[i]["is_ellipse"] = 0
        strokes_topology[i]["old_index"] = i
        strokes_topology[i]["plane_point"] = []
        strokes_topology[i]["plane_normal"] = []

# remove all primitive_type -1
def clean_strokes_topology(strokes_topology, intersections, only_markers=True):
    del_strokes = []
    for i in reversed(range(len(strokes_topology))):
        if only_markers:
            if strokes_topology[i]["primitive_type"] == -2 or \
                isinstance(strokes_topology[i]["points2D"], dict):
                del_strokes.append(strokes_topology[i]["old_index"])
                del strokes_topology[i]
        else:
            if strokes_topology[i]["primitive_type"] == -2 or \
                (strokes_topology[i]["primitive_type"] == 0 and
                 not strokes_topology[i]["depth_assigned"]):
                del_strokes.append(strokes_topology[i]["old_index"])
                del strokes_topology[i]

    del_intersections = []
    for i in sorted(del_strokes, reverse=True):
        for inter_idx, inter in enumerate(intersections):
            if inter["strokes_indices"][0] - 1 == i or \
                    inter["strokes_indices"][1] - 1 == i:
                del_intersections.append(inter_idx)
    for i in sorted(del_intersections, reverse=True):
        del intersections[i]

def inter_is_unique(inter, inter_list):

    line_centroid = inter["straight_linestring"].centroid
    line_idx = inter["straight_line_idx"]

    for snd_inter in inter_list:
        snd_inter_centroid = snd_inter["straight_linestring"].centroid
        snd_line_idx = snd_inter["straight_line_idx"]
        if snd_line_idx == line_idx and \
            np.isclose(line_centroid.distance(snd_inter_centroid), 0.0):
            return False
    return True

def read_stroke_clusters(cluster_file_name, cluster_geometry_file_name,
                         offset_file_name):

    clusters = np.loadtxt(cluster_file_name, dtype=int)
    clusters[:, 0] -= 1
    offset = np.loadtxt(offset_file_name, dtype=float)
    # read cluster_geometry
    nb_clusters = np.max(clusters[:, 1]) + 1
    cluster_lines = []
    for i in range(nb_clusters):
        cluster_lines.append([])
    file = open(cluster_geometry_file_name)
    lines = file.readlines()
    for line_idx in range(len(lines)):
        lines[line_idx] = lines[line_idx].replace("\t", " ")
        lines[line_idx] = lines[line_idx].replace("\n", " ")
    file.close()

    # cluster_scap_ids[i] contains the scap_id of the i-th cluster
    cluster_scap_ids = np.ones(nb_clusters, dtype=int)
    cluster_scap_ids.fill(-1)
    line_idx = 0
    while line_idx < len(lines):
        line = lines[line_idx]
        if line[0] == "{": # new cluster line begins
            # get cluster line index
            line_idx += 1
            line = lines[line_idx]
            line_indices = line.split()
            scap_id = int(line_indices[0][1:]) # without the '#'
            cluster_id = int(line.split()[-1])
            if cluster_id < nb_clusters:
                cluster_scap_ids[cluster_id] = scap_id

            polyline = []
            # next line can be skipped
            line_idx += 2
            line = lines[line_idx]
            while line[0] != "}":
                point = np.array(line.split()[:2], dtype=float)
                polyline.append(point + offset)
                line_idx += 1
                line = lines[line_idx]
            polyline = np.array(polyline)
            if cluster_id < nb_clusters:
                cluster_lines[cluster_id] = polyline
        else:
            line_idx += 1
    # cluster_sets[i] contains the stroke-indices of the i-th cluster
    cluster_sets = []
    for i in range(nb_clusters):
        cluster_sets.append([])
    for c in clusters:
        cluster_sets[c[1]].append(c[0])

    return cluster_sets, cluster_scap_ids.tolist(), cluster_lines

def regroup_strokes(strokes_topology, clusters,
                    cluster_scap_ids, cluster_lines, intersections):

    old_new_stroke_mapping = []
    # assign default value to all strokes
    for s_idx in range(len(strokes_topology)):
        strokes_topology[s_idx]["cluster_id"] = -1
        strokes_topology[s_idx]["scap_id"] = -1
    # assign real clusters
    for cluster_idx, cluster in enumerate(clusters):
        for c in cluster:
            strokes_topology[c]["cluster_id"] = cluster_idx
            strokes_topology[c]["scap_id"] = cluster_scap_ids[cluster_idx]
            strokes_topology[c]["points2D"] = cluster_lines[cluster_idx].tolist()

    # insert polylines at the last curve of the cluster
    del_strokes = []
    for cluster_idx, cluster in enumerate(clusters):
        if len(cluster) < 2:
            continue
        max_cluster_id = np.max(cluster)
        max_pressure = strokes_topology[max_cluster_id]["mean_pressure"]
        for c in cluster:
            if c != max_cluster_id:
                max_pressure = max(max_pressure, strokes_topology[c]["mean_pressure"])
                del_strokes.append(c)
                old_new_stroke_mapping.append([strokes_topology[c]["old_index"], strokes_topology[max_cluster_id]["old_index"]])
                # redirect old intersection indices to new one
                for inter_idx, inter in enumerate(intersections):
                    if inter["strokes_indices"][0] - 1 == strokes_topology[c]["old_index"]:
                        intersections[inter_idx]["strokes_indices"][0] = strokes_topology[max_cluster_id]["old_index"] + 1
                    if inter["strokes_indices"][1] - 1 == strokes_topology[c]["old_index"]:
                        intersections[inter_idx]["strokes_indices"][1] = strokes_topology[max_cluster_id]["old_index"] + 1
        strokes_topology[max_cluster_id]["mean_pressure"] = max_pressure

    for i in sorted(del_strokes, reverse=True):
        del strokes_topology[i]

    return old_new_stroke_mapping

def remove_markers(strokes_topology, intersections):
    strokes_array = convert_strokes_topology_to_strokes_array(strokes_topology)
    del_strokes = []
    old_indices = []
    for s_idx, s in enumerate(strokes_topology):
        if s["primitive_type"] == 1:
            # for curves, the length2D has yet to be established
            points = strokes_array[s_idx]
            bez = [
                np.array(fitCurves.generate_bezier_without_tangents(points))][0]
            # interpolate curve
            interp_points = []
            for t in np.arange(0.0, 1.05, 0.05):
                interp_points.append(bz.q(bez, t))
            curve = LineString(np.array(interp_points))
            s["length2D"] = curve.length
            if s["length2D"] < 5.0*s["accuracy_radius"]:
                del_strokes.append(s_idx)
                old_indices.append(s["old_index"])
        elif len(strokes_array[s_idx]) < 2 or s["accuracy_radius"] is None:
            del_strokes.append(s_idx)
            old_indices.append(s["old_index"])

    for i in sorted(del_strokes, reverse=True):
        del strokes_topology[i]

    del_intersections = []
    for i in sorted(old_indices, reverse=True):
        for inter_idx, inter in enumerate(intersections):
            if inter["strokes_indices"][0] - 1 == i or \
                    inter["strokes_indices"][1] - 1 == i:
                del_intersections.append(inter_idx)
    del_intersections = np.unique(del_intersections)
    for i in sorted(del_intersections, reverse=True):
        del intersections[i]

def check_is_ellipse(points, acc_radius):
    ellipse, u, s = fit_ellipse(points)

    if ellipse is None:
        return None, False

    if np.isclose(s[0], 0.0) or np.isclose(s[1], 0.0):
        return None, False
    # check scaling to detect a line
    if s[0] / s[1] > 10.0 or s[1] / s[0] > 10.0:
        return None, False

    poly = LineString(points).buffer(acc_radius)
    if len(list(poly.interiors)) == 0:
        return None, False

    return ellipse, True


def fit_ellipses(strokes_topology, strokes_array):

    for stroke_idx, stroke in enumerate(strokes_topology):
        strokes_topology[stroke_idx]["is_ellipse"] = 0
        if stroke["primitive_type"] != 1:
            continue
        points = strokes_array[stroke_idx]
        ellipse, is_ellipse = check_is_ellipse(points, stroke["accuracy_radius"])

        if not is_ellipse:
            continue

        strokes_topology[stroke_idx]["is_ellipse"] = is_ellipse
        # fit bezier
        beziers = fitCurves.fitCurve(ellipse, 1.0)
        strokes_topology[stroke_idx]["bezier_cps"] = beziers
        bezier_linestring = []
        for bez in beziers:
            bezier_linestring.extend(bezier_interpolation(bez, 10))
        strokes_topology[stroke_idx]["bezier_linestring"] = LineString(np.array(bezier_linestring))


def fit_ellipse(stroke):

    if len(stroke) == 0:
        return None, None, None
    x = stroke[:, 0].copy()
    y = stroke[:, 1].copy()
    xmean, ymean = x.mean(), y.mean()
    x -= xmean
    y -= ymean
    U, S, V = np.linalg.svd(np.stack((x, y)))
    N = x.shape[0]

    tt = np.linspace(0, 2 * np.pi, 1000)
    circle = np.stack((np.cos(tt), np.sin(tt)))  # unit circle
    transform = np.sqrt(2.0 / float(N)) *np.matmul(U, np.diag(S))  # transformation matrix

    fit = transform.dot(circle) + np.array([[xmean], [ymean]])
    fit = np.transpose(fit)
    return fit, U, S

def get_im_dim(strokes_array):
    x0 = 10000.0
    x1 = -10000.0
    y0 = 10000.0
    y1 = -10000.0
    for s in strokes_array:
        if len(s) > 1:
            x0 = np.minimum(x0, np.min(s[:, 0]))
            x1 = np.maximum(x1, np.max(s[:, 0]))
            y0 = np.minimum(y0, np.min(s[:, 1]))
            y1 = np.maximum(y1, np.max(s[:, 1]))
    return [x0, x1, y0, y1]

def recompute_im_vec(cam_params):

    lambdas = np.ones(1, dtype=float)
    im_dim = cam_params["im_dim"]
    principal_point = np.array([(im_dim[0] + im_dim[1]) / 2.0,
                                (im_dim[0] + im_dim[1]) / 2.0])
    im_center = np.array(
        project_to_3d([principal_point], [lambdas], cam_params)[0])
    cam_pos = cam_params["C"]
    im_vec = im_center - cam_pos
    im_vec /= np.linalg.norm(im_vec)
    cam_params["view_dir"] = im_vec

def project_to_3d(stroke_2d, lambdas, cam_params):

    r_t_inv = cam_params["R_t_inv"]
    k_inv = cam_params["K_inv"]
    stroke_3d_approximated = []
    for p_idx, p in enumerate(stroke_2d):
        u, v = p[0], p[1]

        p_cam = np.dot(k_inv, np.array([[u],[v],[1.0]]))
        p_cam *= lambdas[p_idx]
        p_cam = np.expand_dims(p_cam, 0)

        p_world = np.ones(shape=(4,1))
        p_world[:3] = p_cam
        p_world = np.dot(r_t_inv, p_world)
        p_world[:] /= p_world[3]
        p_world = p_world[:3]
        p_world = np.transpose(p_world)

        stroke_3d_approximated.append(p_world[0])
    return stroke_3d_approximated

def get_stroke_planes(strokes_topology, intersections, strokes_array_3d):
    for stroke_idx, stroke in enumerate(strokes_topology):

        strokes_topology[stroke_idx]["planes"] = []
        strokes_topology[stroke_idx]["confident_curves"] = []
        stroke_3d = strokes_array_3d[stroke_idx]
        if len(stroke_3d) == 0:
            continue

        for inter_idx, inter in enumerate(intersections):
            if inter["collinear"] == 1:
                continue
            if not inter["is_active"] or inter["is_active"] is None:
                continue
            other_stroke_idx = -1
            if inter["strokes_indices"][0] == stroke_idx and \
                    inter["strokes_indices"][1] < stroke_idx:
                other_stroke_idx = inter["strokes_indices"][1]
            if inter["strokes_indices"][1] == stroke_idx and \
                    inter["strokes_indices"][0] < stroke_idx:
                other_stroke_idx = inter["strokes_indices"][0]
            if other_stroke_idx == -1:
                continue
            other_stroke_3d = strokes_array_3d[other_stroke_idx]
            if len(other_stroke_3d) == 0 or \
                    strokes_topology[other_stroke_idx]["primitive_type"] != 0:
                continue
            plane_point = np.array(inter["coordinates3D"])
            vec_1 = stroke_3d[1] - stroke_3d[0]
            if np.isclose(0.0, np.linalg.norm(vec_1)):
                continue
            vec_1 /= np.linalg.norm(vec_1)
            vec_2 = other_stroke_3d[1] - other_stroke_3d[0]
            if np.isclose(0.0, np.linalg.norm(vec_2)):
                continue
            vec_2 /= np.linalg.norm(vec_2)
            if np.abs(np.dot(vec_1, vec_2)) > 0.94:
                continue
            plane_normal = np.cross(vec_1, vec_2)
            plane_normal /= np.linalg.norm(plane_normal)
            strokes_topology[stroke_idx]["planes"].append({"plane_point": plane_point,
                                                           "plane_normal": plane_normal,
                                                           "plane_strokes": [stroke_3d, other_stroke_3d],
                                                           "plane_strokes_indices": [stroke_idx, other_stroke_idx]})

# line is 3d
# p_start, p_end: 2d
def line_substring_points(line, p_start, p_end, cam_params):
    projected_line = np.array(project_3d_stroke(np.array(line), cam_params))
    # go through the projected points and simply get the 2D sub-polyline covered
    # by [p_start, p_end]
    min_start_idx = 0
    min_end_idx = 0
    min_start_dist = p_start.distance(shapely.geometry.Point(projected_line[0]))
    min_end_dist = p_end.distance(shapely.geometry.Point(projected_line[0]))
    for i in range(1, len(projected_line)):
        p = shapely.geometry.Point(projected_line[i])
        start_dist = p.distance(p_start)
        end_dist = p.distance(p_end)
        if start_dist < min_start_dist:
            min_start_dist = start_dist
            min_start_idx = i
        if end_dist < min_end_dist:
            min_end_dist = end_dist
            min_end_idx = i

    # return the corresponding part from line
    seg_length = 0.0
    for i in range(min_start_idx, min_end_idx):
        seg_length += np.linalg.norm(line[i+1]-line[i])
    if len(line[min_start_idx:min_end_idx+1]) < 2:
        min_start_idx = max(0, min_start_idx-1)
        min_end_idx = min(len(line)-1, min_end_idx+1)
    return line[min_start_idx:min_end_idx+1], seg_length

def dist_to_plane(p, plane):
    # http://mathworld.wolfram.com/Point-PlaneDistance.html
    w = p - plane["plane_point"]
    dist_ratio_0 = np.abs(np.dot(plane["plane_normal"], w))
    # check if dist_ratio different at a different point
    w = p - plane["plane_strokes"][0][0]
    dist_ratio_1 = np.abs(np.dot(plane["plane_normal"], w))
    w = p - plane["plane_strokes"][0][1]
    dist_ratio_2 = np.abs(np.dot(plane["plane_normal"], w))
    min_dist = min(dist_ratio_0, min(dist_ratio_1, dist_ratio_2))

    return min_dist

def detect_non_planar_candidates(strokes_topology, include_curves=False,
                                 curve_indices=[]):

    extremity_length = 0.2
    for s_idx, s in enumerate(strokes_topology):
        strokes_topology[s_idx]["non_planar_curve_candidate"] = False
        strokes_topology[s_idx]["non_planar_curve"] = False
        if s["primitive_type"] != 1 or s["is_ellipse"]:
            continue
        if len(curve_indices) > 0 and not s_idx in curve_indices:
            continue
        # check if 2 tangent intersections at the extremities
        curve_start = substring(s["linestring"], 0.0, extremity_length, normalized=True)
        curve_end = substring(s["linestring"], 1.0-extremity_length, 1.0, normalized=True)
        start_covered = False
        end_covered = False
        start_indices = []
        end_indices = []
        for inter_idx, inter in enumerate(s["curve_intersections"]):
            if include_curves:
                if inter["curve_curve_intersection"] and strokes_topology[inter["straight_line_idx"]]["is_ellipse"]:
                    continue
            else:
                if inter["curve_curve_intersection"]:
                    continue

            if inter["is_tangential"] or inter["is_tangential_extended"]:
                if inter["straight_linestring"].intersects(curve_start.buffer(s["accuracy_radius"])):
                    start_covered = True
                    start_indices.append(inter_idx)
                elif inter["straight_linestring"].intersects(curve_end.buffer(s["accuracy_radius"])):
                    # elif because we do not want curve_start and curve_end be
                    # covered by the same curve
                    end_covered = True
                    end_indices.append(inter_idx)
        if start_covered and end_covered:
            strokes_topology[s_idx]["non_planar_curve_candidate"] = True
            strokes_topology[s_idx]["non_planar_start_indices_candidates"] = start_indices
            strokes_topology[s_idx]["non_planar_end_indices_candidates"] = end_indices
            # todo: think about what to do with multiple start_indices and end_indices
            # solution: check for any combination of start_indices, end_indices
            # if a planar solution can be formed. If that's the case, then don't
            # consider this curve being a planar curve

            non_planar_start_indices = []
            non_planar_end_indices = []
            planar_solution_found = False
            for start_idx in start_indices:
                for end_idx in end_indices:
                    # form plane with the 2 tangents and the start_point, then check if
                    # end_point is in the plane
                    start_inter = s["curve_intersections"][start_idx]
                    end_inter = s["curve_intersections"][end_idx]
                    if include_curves:
                        if strokes_topology[start_inter["straight_line_idx"]]["primitive_type"] == 1\
                            or strokes_topology[end_inter["straight_line_idx"]]["primitive_type"] == 1:
                            non_planar_start_indices.append(start_idx)
                            non_planar_end_indices.append(end_idx)
                            continue

                    start_tangent = np.array(start_inter["positional_constraint_3d"])
                    start_point = start_tangent[0]
                    start_tangent = start_tangent[-1] - start_tangent[0]
                    end_tangent = np.array(end_inter["positional_constraint_3d"])
                    end_point = end_tangent[0]
                    end_tangent = end_tangent[-1] - end_tangent[0]
                    if np.isclose(np.linalg.norm(start_tangent), 0.0):
                        continue
                    if np.isclose(np.linalg.norm(end_tangent), 0.0):
                        continue
                    start_tangent /= np.linalg.norm(start_tangent)
                    end_tangent /= np.linalg.norm(end_tangent)
                    # check if both tangents are too colinear
                    if np.abs(np.dot(start_tangent, end_tangent)) > 0.94:
                        planar_solution_found = True # todo: correct
                        continue
                    plane_normal = np.cross(start_tangent, end_tangent)
                    plane_normal /= np.linalg.norm(plane_normal)
                    plane = {"plane_point": start_point,
                            "plane_normal": plane_normal,
                            "plane_strokes": [np.array([start_point, start_point+start_tangent]),
                                              np.array([start_point, start_point+end_tangent])],
                             "plane_strokes_indices": [0, 0]}

                    dist_ratio = dist_to_plane(end_point, plane)

                    if dist_ratio > 0.010:
                        non_planar_start_indices.append(start_idx)
                        non_planar_end_indices.append(end_idx)
                    else:
                        planar_solution_found = True
            strokes_topology[s_idx]["non_planar_start_indices"] = non_planar_start_indices
            strokes_topology[s_idx]["non_planar_end_indices"] = non_planar_end_indices
            strokes_topology[s_idx]["non_planar_curve"] = not planar_solution_found

def reject_intersection(inter, strokes_topology, reject=True):
    is_ellipse = strokes_topology[inter["curved_line_idx"]]["is_ellipse"]
    other_is_ellipse = strokes_topology[inter["straight_line_idx"]]["is_ellipse"]
    if inter["curve_curve_intersection"]:
        if (is_ellipse and other_is_ellipse) and not (inter["is_tangential"] or inter["is_tangential_extended"]):
            return True
        if is_ellipse and not other_is_ellipse:
            return False
        if other_is_ellipse:
            return False
        if inter["is_tangential"] or inter["is_tangential_extended"]:

            t_params = inter["straight_t_param"]
            first_overlap = abs(max(t_params) - min(t_params))
            t_params = inter["curved_t_param"]
            snd_overlap = abs(max(t_params) - min(t_params))
            if reject:
                if first_overlap > 0.3 and snd_overlap > 0.3:
                    return True
    return False

def filter_out_curve_intersections(curve_intersections, strokes_topology):

    del_indices = []
    for inter_idx, inter in enumerate(curve_intersections):
        if reject_intersection(inter, strokes_topology):
            del_indices.append(inter_idx)

    for del_idx in sorted(del_indices, reverse=True):
        del curve_intersections[del_idx]

# returns two lists:
# parent_list
# child_list
# the two lists are of the same size
# children are unique in child_list, parents can be parents to more than one child
def get_absorbed_strokes(strokes_topology):
    parent_list = []
    child_list = []
    for s_idx, s in enumerate(strokes_topology):
        if s["primitive_type"] != 1:
            continue

        # get curves absorbed by multiple lines
        tan_intersections = []
        tan_overlaps = []
        tan_s_indices = []
        overlap_interval = interval()
        for inter in s["curve_intersections"]:
            if inter["curve_curve_intersection"] and \
                    (inter["is_tangential"] or inter["is_tangential_extended"]):
                tan_intersections.append(inter)
                tan_overlaps.append(inter["curved_t_param"])
                overlap_interval = overlap_interval | interval[min(inter["curved_t_param"]), max(inter["curved_t_param"])]
                tan_s_indices.append(inter["straight_line_idx"])
        if len(overlap_interval) == 0:
            continue
        overlap = 0.0
        for tmp_inter in overlap_interval.components:
            overlap += tmp_inter[0][1] - tmp_inter[0][0]
        if overlap < 0.9:
            continue
        # parent_idx = curve with the biggest overlap
        separate_overlaps = []
        for over in tan_overlaps:
            separate_overlaps.append(abs(max(over)-min(over)))
        parent_idx = tan_s_indices[np.argmax(separate_overlaps)]
        if parent_idx in child_list:
            parent_idx_idx = np.argwhere(np.array(child_list)[:] == parent_idx)[0][0]
            parent_idx = parent_list[parent_idx_idx]
        parent_list.append(parent_idx)
        child_list.append(s_idx)
    return parent_list, child_list

def cluster_intersecting_strokes_networkx_v2(stroke_indices, strokes_topology):

    def first_elem(list):
        return list[0]

    # build adjacency matrix
    adj_mat = np.zeros(shape=[len(stroke_indices), len(stroke_indices)],
                       dtype=int)
    for s_idx_1_idx, s_idx_1 in enumerate(stroke_indices):
        intersecting_s_ids = [inter["straight_line_idx"] for inter in
                              strokes_topology[s_idx_1]["curve_intersections"]]
        for inter in strokes_topology[s_idx_1]["curve_intersections"]:
            s_idx_2 = np.argwhere(np.array(stroke_indices) == inter["straight_line_idx"]).flatten()
            if len(s_idx_2) > 0:
                s_idx_2 = s_idx_2[0]
                adj_mat[s_idx_1_idx][s_idx_2] = 1
                adj_mat[s_idx_2][s_idx_1_idx] = 1

    # enforce symmetry
    for i in range(len(stroke_indices)):
        for j in range(len(stroke_indices)):
            if adj_mat[i][j] == 1:
                adj_mat[j][i] = 1
    adj_graph = nx.from_numpy_matrix(adj_mat)
    clusters = [sorted(stroke_indices[list(c)].tolist())
                for c in sorted(nx.connected_components(adj_graph), key=len, reverse=True)]
    clusters = sorted(clusters, key=first_elem)
    return clusters
