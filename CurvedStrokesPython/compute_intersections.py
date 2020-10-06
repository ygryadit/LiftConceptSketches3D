import tools
import numpy as np
import shapely.geometry

def compute_curve_2D_intersections_single_curve(curve_idx,
                                                strokes_array, strokes_array_3d,
                                                lines_array, lines_array_3d,
                                                strokes_topology, intersections,
                                                cc_intersections,
                                                tangent_sensibility=0.4,
                                                VERBOSE=0, ellipse_clusters=[],
                                                reject=True):
    curve_intersections = []
    stroke_start_intersections = []
    stroke_end_intersections = []
    stroke_start_covered = False
    stroke_end_covered = False
    curved_line = lines_array[curve_idx]
    acc_radius = strokes_topology[curve_idx]["accuracy_radius"]
    # intersections with previous strokes
    for prev_stroke_idx, prev_stroke in enumerate(
            strokes_array[:curve_idx]):
        # but only with straight lines
        if not(strokes_topology[prev_stroke_idx]["primitive_type"] == 0 or \
               strokes_topology[prev_stroke_idx]["primitive_type"] == 1):
            continue

        if cc_intersections:
            if strokes_topology[prev_stroke_idx]["primitive_type"] == 0 and \
                    not strokes_topology[prev_stroke_idx]["depth_assigned"]:
                continue
        else:
            if strokes_topology[prev_stroke_idx]["primitive_type"] != 0 or \
                    not strokes_topology[prev_stroke_idx]["depth_assigned"]:
                continue
        curve_curve_intersection = \
            strokes_topology[prev_stroke_idx]["primitive_type"] == 1

        prev_line = lines_array[prev_stroke_idx]
        if curved_line.buffer(acc_radius).intersects(prev_line):
            # actual intersections - on the straight line
            intersection_point_list = tools.get_intersection_list(
                curved_line,
                prev_line, acc_radius)

            # get corresponding parts from the current curved stroke
            for intersected_line_idx, intersected_line in enumerate(
                    intersection_point_list):
                tmp_line_straight = shapely.geometry.LineString(
                    intersected_line)
                # intersected line segments on the curved line
                prior_intersection_point_list = tools.get_intersection_list(
                    tmp_line_straight,
                    curved_line, acc_radius)
                for p_idx, p in enumerate(prior_intersection_point_list):
                    # check if the two lines are tangential
                    if len(p) == 0:
                        continue
                    tmp_line_curved = shapely.geometry.LineString(p)
                    is_tangential = tools.is_tangential(
                        tmp_line_curved, tmp_line_straight,
                        tangent_sensibility)
                    is_tangential_extended = tools.is_tangential_extended(
                        tmp_line_curved, tmp_line_straight,
                        tangent_sensibility)
                    is_trihedral_end = tools.is_trihedral_end(curve_idx,
                                                              tmp_line_curved,
                                                              tmp_line_straight,
                                                              intersections,
                                                              acc_radius,
                                                              strokes_array,
                                                              strokes_topology)
                    # also take into account stroke extremities
                    inter_start_covered = tools.stroke_start_covered(
                        curve_idx,
                        tmp_line_curved, strokes_array, acc_radius)
                    inter_end_covered = tools.stroke_end_covered(curve_idx,
                                                                 tmp_line_curved,
                                                                 strokes_array,
                                                                 acc_radius)
                    inter_dict = {}
                    if is_tangential or is_tangential_extended or \
                            is_trihedral_end or inter_start_covered or \
                            inter_end_covered or curve_curve_intersection:
                        # get intersected line in 3d
                        prev_line_3d = lines_array_3d[prev_stroke_idx]
                        if not curve_curve_intersection and prev_line_3d.is_empty:
                            continue
                        # filter out non-tangential intersections in the middle of
                        # the curve, but only with straight lines
                        if not strokes_topology[curve_idx]["is_ellipse"] and \
                                not curve_curve_intersection and \
                            not (inter_start_covered or inter_end_covered) and \
                            (not (is_tangential or is_tangential_extended)):
                            continue

                        # get t_params for curved line and for straight line
                        curve_t_params = []
                        curve_t_params.append(tools.get_closest_param_beziers(
                            strokes_topology[curve_idx]["bezier_cps"],
                            np.array(tmp_line_straight.interpolate(0.0, normalized=True))))
                        curve_t_params.append(tools.get_closest_param_beziers(
                            strokes_topology[curve_idx]["bezier_cps"],
                            np.array(tmp_line_straight.centroid)))
                        curve_t_params.append(tools.get_closest_param_beziers(
                            strokes_topology[curve_idx]["bezier_cps"],
                            np.array(tmp_line_straight.interpolate(1.0, normalized=True))))
                        straight_t_params = []
                        if curve_curve_intersection:
                            straight_t_params.append(
                                tools.get_closest_param_beziers(
                                    strokes_topology[prev_stroke_idx]["bezier_cps"],
                                    np.array(tmp_line_curved.interpolate(0.0,
                                                                           normalized=True))))
                            straight_t_params.append(tools.get_closest_param_beziers(
                                strokes_topology[prev_stroke_idx]["bezier_cps"],
                                np.array(tmp_line_curved.centroid)))
                            straight_t_params.append(
                                tools.get_closest_param_beziers(
                                    strokes_topology[prev_stroke_idx]["bezier_cps"],
                                    np.array(tmp_line_curved.interpolate(1.0,
                                                                         normalized=True))))

                        inter_dict = {"curved_line_idx": curve_idx,
                                      "straight_line_idx": prev_stroke_idx,
                                      "curved_linestring": tmp_line_curved,
                                      "straight_linestring": tmp_line_straight,
                                      "curved_t_param": curve_t_params,
                                      "straight_t_param": straight_t_params,
                                      "is_tangential": is_tangential,
                                      "is_tangential_extended": is_tangential_extended,
                                      "is_trihedral": is_trihedral_end,
                                      "accuracy_radius":
                                          strokes_topology[curve_idx][
                                              "accuracy_radius"],
                                      "curve_curve_intersection":
                                          curve_curve_intersection}
                        if tools.reject_intersection(inter_dict, strokes_topology,
                                                     reject=reject):
                            continue
                        stroke_start_covered = stroke_start_covered or \
                                               ((
                                                        is_tangential or is_tangential_extended or is_trihedral_end)
                                                and inter_start_covered)
                        stroke_end_covered = stroke_end_covered or \
                                             ((
                                                      is_tangential or is_tangential_extended or is_trihedral_end)
                                              and inter_end_covered)
                    if is_tangential or is_tangential_extended or \
                            is_trihedral_end or curve_curve_intersection:
                        if tools.inter_is_unique(inter_dict,
                                                 curve_intersections):
                            curve_intersections.append(inter_dict)
                    elif inter_start_covered:
                        if tools.inter_is_unique(inter_dict,
                                                 stroke_start_intersections):
                            stroke_start_intersections.append(inter_dict)
                    elif inter_end_covered:
                        if tools.inter_is_unique(inter_dict,
                                                 stroke_end_intersections):
                            stroke_end_intersections.append(inter_dict)
    if not stroke_start_covered:
        if len(stroke_start_intersections) > 0:
            curve_intersections.extend(stroke_start_intersections)
    if not stroke_end_covered:
        if len(stroke_end_intersections) > 0:
            curve_intersections.extend(stroke_end_intersections)
    return curve_intersections

# function specifically designed to reject non-planar clustered curves
def compute_curve_2D_intersections_single_curve_input_linestring_simplified(curve_idx,
                                                                 curved_line,
                                                                 acc_radius,
                                                                 strokes_array,
                                                                 lines_array, lines_array_3d,
                                                                 strokes_topology,
                                                                 tangent_sensibility=0.4):
    curve_intersections = []
    curved_line = curved_line
    # intersections with previous strokes
    for prev_stroke_idx, prev_stroke in enumerate(
            strokes_array[:curve_idx]):
        # but only with straight lines
        if strokes_topology[prev_stroke_idx]["primitive_type"] != 0 or \
                not strokes_topology[prev_stroke_idx]["depth_assigned"]:
            continue

        prev_line = lines_array[prev_stroke_idx]
        if curved_line.buffer(acc_radius).intersects(prev_line):
            # actual intersections - on the straight line
            intersection_point_list = tools.get_intersection_list(
                curved_line,
                prev_line, acc_radius)

            # get corresponding parts from the current curved stroke
            for intersected_line_idx, intersected_line in enumerate(
                    intersection_point_list):
                tmp_line_straight = shapely.geometry.LineString(
                    intersected_line)
                # intersected line segments on the curved line
                prior_intersection_point_list = tools.get_intersection_list(
                    tmp_line_straight,
                    curved_line, acc_radius)
                for p_idx, p in enumerate(prior_intersection_point_list):
                    # check if the two lines are tangential
                    if len(p) == 0:
                        continue
                    tmp_line_curved = shapely.geometry.LineString(p)
                    is_tangential = tools.is_tangential(
                        tmp_line_curved, tmp_line_straight,
                        tangent_sensibility)
                    is_tangential_extended = tools.is_tangential_extended(
                        tmp_line_curved, tmp_line_straight,
                        tangent_sensibility)

                    inter_dict = {}
                    if is_tangential or is_tangential_extended:
                        # get intersected line in 3d
                        prev_line_3d = lines_array_3d[prev_stroke_idx]
                        if prev_line_3d.is_empty:
                            continue

                        inter_dict = {"straight_line_idx": prev_stroke_idx,
                                      "straight_linestring": tmp_line_straight,
                                      "is_tangential": is_tangential,
                                      "is_tangential_extended": is_tangential_extended,
                                      "curve_curve_intersection": False}

                    if is_tangential or is_tangential_extended:
                        if tools.inter_is_unique(inter_dict,
                                                 curve_intersections):
                            curve_intersections.append(inter_dict)
    return curve_intersections
