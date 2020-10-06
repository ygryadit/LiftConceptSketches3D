import numpy as np
from math import acos, sqrt, pi, exp
import polyscope as ps

ALPHA = np.deg2rad(15.0)
SIGMA_1 = (1.0 - np.cos(ALPHA))/3.0
SIGMA_2 = (np.cos(np.deg2rad(90 - 15)))/3.0

# a curve has a geometry, a plane_point and a plane_normal
class Curve3D:

    def __init__(self, geometry=None, plane_point=None, plane_normal=None):
        self.geometry = geometry
        self.plane_point = plane_point
        self.plane_normal = plane_normal

def line_plane_collision(plane_normal, plane_point, ray_dir, ray_p):

    ndotu = plane_normal.dot(ray_dir)
    if np.isclose(abs(ndotu), 0.0):
        return None

    w = ray_p - plane_point
    si = -plane_normal.dot(w) / ndotu
    return w + si * ray_dir + plane_point

# v1 and v2 should be normalized
# returns closest points on line 1 and on line 2
def line_line_collision(p1, v1, p2, v2):
    v3 = np.cross(v1, v2)
    v3 /= np.linalg.norm(v3)

    rhs = p2 - p1
    lhs = np.array([v1, -v2, v3]).T

    t_solutions = np.linalg.lstsq(lhs, rhs, rcond=None)
    t1 = t_solutions[0][0]
    t2 = t_solutions[0][1]

    closest_line_1 = p1 + t1*v1
    closest_line_2 = p2 + t2*v2
    return [closest_line_1, closest_line_2]

def get_ellipse_eccentricity(ellipse_3d, plane_u, plane_v):

    if np.dot(plane_u, plane_v) < 0:
        plane_v *= -1.0
    projected_ellipse = np.array([[np.dot(plane_u, p), np.dot(plane_v, p)] for p in ellipse_3d])
    x = projected_ellipse[:, 0]
    y = projected_ellipse[:, 1]
    xmean, ymean = x.mean(), y.mean()
    x -= xmean
    y -= ymean
    _, scale, _ = np.linalg.svd(np.stack((x, y)))

    major_axis_length = np.max(scale)
    minor_axis_length = np.min(scale)
    a = major_axis_length / 2
    b = minor_axis_length / 2
    ecc = np.sqrt(np.square(a) - np.square(b)) / a
    return ecc

def distance_point_to_line(p, line_p, line_v):
    dist = np.linalg.norm(np.cross(line_p - p, line_v))/np.linalg.norm(line_v)
    return dist

# 3D polyline
def distance_point_to_polyline(p, polyline):

    print("distance normal")
    p = np.array(p)
    min_dist = 1000.0
    for p_id in range(len(polyline)-1):
        # check that p is within [a, b]
        a = np.array(polyline[p_id])
        b = np.array(polyline[p_id+1])
        ab = b - a
        ab /= np.linalg.norm(ab)
        ap = p - a
        ap /= np.linalg.norm(ap)
        bp = p - b
        bp /= np.linalg.norm(bp)

        if np.dot(ap, ab) <= 0.0:
            min_dist = min(min_dist, np.linalg.norm(p-a))
        elif np.dot(bp, ab) >= 0.0:
            min_dist = min(min_dist, np.linalg.norm(p-b))
        else:
            min_dist = min(min_dist, distance_point_to_line(p, a, ab))

    return min_dist


def lineseg_dist(p, a, b):

    p = np.array([p])
    p = np.repeat(p, len(a), axis=0)
    # normalized tangent vector
    if np.any(np.isclose(0.0, np.linalg.norm(b - a, axis=-1).reshape(-1, 1))):
        return []
    d = np.divide(b - a, np.linalg.norm(b - a, axis=-1).reshape(-1, 1))

    # signed parallel distance components
    s = np.einsum("ij, ij->i", a-p, d)
    t = np.einsum("ij, ij->i", p-b, d)

    # clamped parallel distance
    h = np.maximum.reduce([s, t, np.zeros(len(p))])
    c = np.cross(p - a, d)

    return np.hypot(h, np.linalg.norm(c, axis=-1))

def distance_point_to_polyline_vectorized(p, polyline):
    """Cartesian distance from point to line segment

    Edited to support arguments as series, from:
    https://stackoverflow.com/a/54442561/11208892

    Args:
        - p: np.array of single point, shape (2,) or 2D array, shape (x, 2)
        - a: np.array of shape (x, 2)
        - b: np.array of shape (x, 2)
    """
    p = np.array(p)
    polyline = np.array(polyline)
    a = polyline[:polyline.shape[0]-1]
    b = polyline[1:]
    distances = lineseg_dist(p, a, b)
    if len(distances) == 0:
        return -1.0
    return np.min(distances)

def distance_point_to_plane(p, plane_p, plane_n):
    p = np.array(p)
    plane_p = np.array(plane_p)
    plane_n = np.array(plane_n)
    plane_n /= np.linalg.norm(plane_n)
    w = p - plane_p
    dist = np.abs(np.dot(w, plane_n))/np.linalg.norm(plane_n)
    return dist

def angle_line_line(l1, l2):
    vec_1 = l1[-1] - l1[0]
    vec_2 = l2[-1] - l2[0]
    norm_vec_1 = vec_1/np.linalg.norm(vec_1)
    norm_vec_2 = vec_2/np.linalg.norm(vec_2)
    if np.isclose(np.dot(norm_vec_1, norm_vec_2), 1.0):
        return 0.0
    return acos(np.abs(np.dot(norm_vec_1, norm_vec_2)))

def line_3d_length(polyline):
    l = 0.0
    for i in range(len(polyline)-1):
        l += np.linalg.norm(polyline[i+1] - polyline[i])
    return l

def merge_two_line_segments(l1, l2):
    # just merge the closest endpoints
    if np.linalg.norm(l1[0] - l2[0]) < np.linalg.norm(l1[0] - l2[1]):
        return np.array([(l1[0]+l2[0])/2.0, (l1[1]+l2[1])/2.0])
    return np.array([(l1[0]+l2[1])/2.0, (l1[1]+l2[0])/2.0])

def merge_n_line_segments(lines):
    l1 = lines[0]
    for l in lines[1:]:
        l1 = merge_two_line_segments(l1, l)
    return l1

def merge_two_curves(c1, c2):
    merged_geom = [(c1.geometry[p_id]+c2.geometry[p_id])/2.0 for p_id in range(len(c1.geometry))]
    merged_geom = np.array(merged_geom)
    merged_plane_point = (c1.plane_point + c2.plane_point)/2.0

    c1_plane_normal = np.array(c1.plane_normal)/np.linalg.norm(c1.plane_normal)
    c2_plane_normal = np.array(c2.plane_normal)/np.linalg.norm(c2.plane_normal)
    normal_angle = np.dot(c1_plane_normal, c2_plane_normal)
    if normal_angle <= 0.0:
        c2_plane_normal *= -1.0
    if np.isclose(1.0, abs(np.dot(c1_plane_normal, c2_plane_normal))):
        merged_plane_normal = c1_plane_normal
    else:
        merged_plane_normal = (c1_plane_normal + c2_plane_normal)/2.0
        merged_plane_normal /= np.linalg.norm(merged_plane_normal)
    return Curve3D(geometry=merged_geom, plane_point=merged_plane_point,
                   plane_normal=merged_plane_normal)

def merge_n_curves(curves, VERBOSE=False, intersections=[]):
    c1 = curves[0]
    for c in curves[1:]:
        c1 = merge_two_curves(c1, c)
    return c1

def line_cylinder_collisions(line_p, line_v, cylinder_o, cylinder_dir,
                                cylinder_radius):
    # test intersection against infinite body
    line_v /= np.linalg.norm(line_v)
    cylinder_dir /= np.linalg.norm(cylinder_dir)
    d = np.cross(line_v, cylinder_dir)
    e = np.cross(line_p - cylinder_o, cylinder_dir)
    # a*t^2 + b*t + c = 0
    a = np.dot(d, d)
    b = 2.0*np.dot(d, e)
    c = np.dot(e, e) - cylinder_radius**2

    det = b**2 - 4.0*a*c
    if det < 0:
        return []
    if np.isclose(det, 0.0):
        t = -b/(2.0*a)
        return [line_p+t*line_v]
    t_1 = (-b + np.sqrt(det))/(2.0*a)
    t_2 = (-b - np.sqrt(det))/(2.0*a)
    return [line_p + t_1 * line_v, line_p + t_2 * line_v]

def gaussian(x, sigma):
    return exp(-0.5*(x**2/sigma**2))

def compute_axis_alignment(l, axis_label):

    # add axis-alignment score to candidate line
    major_axis = np.zeros(3)
    major_axis[axis_label] = 1.0
    line_vec = l[-1] - l[0]
    line_vec /= np.linalg.norm(line_vec)
    axis_alignment = gaussian(1.0 - abs(np.dot(line_vec, major_axis)), SIGMA_1)
    return axis_alignment

def compute_dot_prod_between_lines(l1, l2):
    l1_vec = l1[-1] - l1[0]
    l1_vec /= np.linalg.norm(l1_vec)
    l2_vec = l2[-1] - l2[0]
    l2_vec /= np.linalg.norm(l2_vec)
    return abs(np.dot(l1_vec, l2_vec))

def compute_gaussian_between_lines(l1, l2):
    l1_vec = l1[-1] - l1[0]
    l1_vec /= np.linalg.norm(l1_vec)
    l2_vec = l2[-1] - l2[0]
    l2_vec /= np.linalg.norm(l2_vec)
    return gaussian(abs(np.dot(l1_vec, l2_vec)), SIGMA_2)

def compute_gaussian_between_line_vec(l1, vec):
    l1_vec = l1[-1] - l1[0]
    l1_vec /= np.linalg.norm(l1_vec)
    vec /= np.linalg.norm(vec)
    return gaussian(abs(np.dot(l1_vec, vec)), SIGMA_2)

def compute_gaussian_between_one_minus_lines(l1, l2):
    l1_vec = l1[-1] - l1[0]
    l1_vec /= np.linalg.norm(l1_vec)
    l2_vec = l2[-1] - l2[0]
    l2_vec /= np.linalg.norm(l2_vec)
    return gaussian(1.0 - abs(np.dot(l1_vec, l2_vec)), SIGMA_2)

def compute_planarity_score(l1, plane_l1, plane_l2):
    plane_vec_1 = plane_l1[-1] - plane_l1[0]
    plane_vec_1 /= np.linalg.norm(plane_vec_1)
    plane_vec_2 = plane_l2[-1] - plane_l2[0]
    plane_vec_2 /= np.linalg.norm(plane_vec_2)
    if abs(np.dot(plane_vec_1, plane_vec_2)) > np.cos(np.deg2rad(80)):
        return 0.0
    plane_normal = np.cross(plane_vec_1, plane_vec_2)
    plane_normal /= np.linalg.norm(plane_normal)
    l1_vec = l1[-1] - l1[0]
    l1_vec /= np.linalg.norm(l1_vec)
    return gaussian(abs(np.dot(l1_vec, plane_normal)), SIGMA_2)

def get_rotation_mat_z(angle):

    return np.array([[np.cos(angle), -np.sin(angle), 0],
                     [np.sin(angle), np.cos(angle), 0],
                     [0, 0, 1]])

def distance_to_bbox(p, bbox):
    p = np.array(p)
    bbox_center = (bbox[:3]+bbox[3:])/2.0
    bbox_dims = bbox[3:] - bbox[:3]
    dx = np.maximum(abs(p[0] - bbox_center[0]) - bbox_dims[0], 0)
    dy = np.maximum(abs(p[1] - bbox_center[1]) - bbox_dims[1], 0)
    dz = np.maximum(abs(p[2] - bbox_center[2]) - bbox_dims[2], 0)
    return np.sqrt(dx*dx + dy*dy + dz*dz)
