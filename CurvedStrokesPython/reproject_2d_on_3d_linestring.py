import tools
import numpy as np
import shapely.geometry

def iter_lin_search(p_2d, line_3d, cam_params):

    line_3d_coords = np.array(line_3d)
    stepsize = 0.001
    t_potentials = np.arange(0.0, 1.0+stepsize, stepsize)
    probes = np.zeros([t_potentials.shape[0], 3])
    probes[:] = line_3d_coords[-1] - line_3d_coords[0]
    t_potentials = np.reshape(t_potentials, [t_potentials.shape[0], 1])
    probes *= t_potentials
    probes[:] += line_3d_coords[0]

    projected_probes = np.array(tools.project_3d_stroke(probes, cam_params))

    distances = np.linalg.norm(projected_probes[:] - p_2d, axis=1)

    return t_potentials[np.argmin(distances)]

def get_3d_constraints_single_curve(curve_intersections,
                           strokes_array,
                           strokes_array_3d,
                           strokes_topology,
                           cam_params,
                           curves_3d=[],
                           cc_intersections=True,
                           only_cc_intersections=False,
                           VERBOSE=0):

    keep_intersections = []
    for inter_idx, inter in enumerate(curve_intersections):
        if not cc_intersections and inter["curve_curve_intersection"]:
            keep_intersections.append(inter_idx)
            continue
        if only_cc_intersections and not inter["curve_curve_intersection"]:
            keep_intersections.append(inter_idx)
            continue
        line_2d = inter["straight_linestring"]
        # if it intersects with a curve, look for the 3D curve in curves_3d
        line_3d = []
        if only_cc_intersections and inter["curve_curve_intersection"]:
            for c_3d in curves_3d:
                if c_3d["curve_idx"] == inter["straight_line_idx"]:
                    line_3d = c_3d["curve_3d"]
        else:
            inter["full_line_3d"] = []
            if len(strokes_array_3d[inter["straight_line_idx"]]) == 0:
                continue
            line_3d = shapely.geometry.LineString(strokes_array_3d[inter["straight_line_idx"]])
            inter["full_line_3d"] = line_3d

        start_point_2d = line_2d.interpolate(0.0, normalized=True)
        end_point_2d = line_2d.interpolate(1.0, normalized=True)

        t_start = iter_lin_search(start_point_2d, line_3d, cam_params)
        t_end = iter_lin_search(end_point_2d, line_3d, cam_params)

        p0 = np.array(tools.line_interpolate(line_3d, t_start))
        p1 = np.array(tools.line_interpolate(line_3d, t_end))

        line_segment_3d = shapely.geometry.LineString([p0, p1])

        if inter["curve_curve_intersection"]:
            inter["curve_3d"] = line_3d
            line_segment_3d, new_length = tools.line_substring_points(line_3d,
                start_point_2d, end_point_2d, cam_params)

        if len(np.array(line_segment_3d)) > 0:
            # remove zero-length intersections
            keep_intersections.append(inter_idx)

        inter["straight_linestring_3d"] = line_segment_3d
        # we want the positional constraint to be as close as possible to the
        # curve, to minmize projection error. Thus, we do not choose the 3D
        # coordinates of the "real" intersection point (matlab)
        inter["positional_constraint_3d"] = line_segment_3d

    return np.array(curve_intersections)[keep_intersections]
