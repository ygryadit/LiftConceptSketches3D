import numpy as np
from copy import deepcopy
import polyscope as ps
import tools_3d
from more_itertools import distinct_combinations

# organized as follows:

# DependencyGraph represents on a high-level the 2D depencies between strokes,
# given by 2D intersections

# Each DependencyNode contains a list of CandidateNodes, which is where the
# different 3D candidate lines are stocked

# Each CandidateNode has a candidate_number, and a list of configurations,
# which represent the different combinations of incoming intersections

class Context3D:
    def __init__(self, sketch_versions=None, dep_node_ids=None, cand_node_ids=None):
        self.sketch_versions = sketch_versions
        self.dep_node_ids = dep_node_ids
        self.cand_node_ids = cand_node_ids

class Intersection:
    def __init__(self, inter_id=None, coords_3d=None, coords_2d=None,
                 tangents_2d=None, tangents_3d=None,
                 stroke_ids=None, inter_params=None, mid_inter_param=None,
                 acc_radius=None, adjacent_inter_ids=None):
        # the inter_id corresponds to the same inter_id used in sketch.intersection_graph
        self.inter_id = inter_id
        self.coords_3d = coords_3d
        self.coords_2d = coords_2d
        self.tangents_2d = tangents_2d
        self.tangents_3d = tangents_3d
        self.tangents_angle_3d = -1.0
        # inter_params contains the mid_inter_params for each stroke
        self.inter_params = inter_params
        self.mid_inter_param = mid_inter_param
        self.stroke_ids = stroke_ids
        self.acc_radius = acc_radius
        self.adjacent_inter_ids = adjacent_inter_ids

class ScoreContainer:
    def __init__(self, line_coverage=-1.0, axis_alignment=-1.0,
                 orthogonality=-1.0, planarity=-1.0, tangentiality=-1.0,
                 foreshortening=-1.0, curve_geom=-1.0, circularity=-1.0):
        self.line_coverage = line_coverage
        self.axis_alignment = axis_alignment
        self.orthogonality = orthogonality
        self.planarity = planarity
        self.tangentiality = tangentiality
        self.foreshortening = foreshortening
        self.curve_geom = curve_geom
        self.circularity = circularity
        self.total_score = 0.0

class CandidateNode:
    def __init__(self, candidate_nb=-1):
        self.candidate_nb = candidate_nb # the n-th candidate of the stroke
        self.geometry = []
        self.plane_point = None
        self.plane_normal = None
        self.intersections = []
        # previous candidates is a list of tuples (dep_node_id, candidate_id)
        self.previous_candidates = []
        # a stroke candidate represents the 3D geometry of a stroke which can
        # belong to several different sketch-versions
        self.sketch_versions = []

        # score function related attributes
        self.axis_alignment = -1.0

        # curve-related information
        self.curve_plane_point = None
        self.curve_plane_normal = None

    # compatibility between self.previous_candidates and cand.previous_candidates
    # compatible if for the same dep_node_id, they have the same cand_id
    def is_compatible(self, cand):
        if len(self.previous_candidates) == 0 or len(cand.previous_candidates) == 0:
            return True
        previous_candidates_2 = np.array(cand.previous_candidates)
        for prev_cand_1 in self.previous_candidates:
            prev_cand_1_loc = np.argwhere(previous_candidates_2[:, 0] == prev_cand_1).flatten()
            if len(prev_cand_1_loc) > 0:
                if previous_candidates_2[prev_cand_1_loc[0]] != prev_cand_1[1]:
                    return False
        return True

    def __str__(self):
        output = "\nCandidateNode: "+str(self.candidate_nb)
        output += "\nsketch_versions: "+str(self.sketch_versions)
        return output

class DependencyNode:
    def __init__(self, stroke_id=-1, in_edges=[], out_edges=[], axis_label=-1,
                 is_curve=False, is_ellipse=False):
        self.stroke_id = stroke_id
        self.in_strokes = in_edges
        self.out_strokes = out_edges

        self.is_assigned = False

        # 3D information about candidate lines
        # the key to a candidate_node is its candidate_nb
        self.candidate_nodes = dict()

        self.axis_label = axis_label
        self.is_curve = is_curve
        self.is_ellipse = is_ellipse

        self.length_2d = 0.0

    def insert_candidate_node(self, cand_node):
        self.candidate_nodes[cand_node.candidate_nb] = cand_node

    def __str__(self):
        output = "\nDependencyNode: "
        output += "\nstroke_id: " + str(self.stroke_id)
        output += "\nin_strokes: " + str(self.in_strokes)
        output += "\nout_strokes: " + str(self.out_strokes)
        for cand_nb in self.candidate_nodes.keys():
            output += str(self.candidate_nodes[cand_nb])
        return output

class DependencyGraph:
    def __init__(self):
        self.dependency_nodes = []
        self.stroke_id_to_dep_node_id = []
        self.sketch_version_counter = 0
        # each entry of sketch_versions_reference corresponds to the description
        # of a sketch_version.
        # This description is a list of triplets.
        # The i-th triplet represents the i-th stroke (in dep_node order, not necessarily original drawing order!)
        # Each triplet contains: (cand_id, [dep_inter_set], cand_score)
        self.sketch_versions_reference = dict()
        self.sketch_versions_scores = dict()
        self.median_length_3d = 0.0
        self.stroke_lengths_3d = dict()
        self.total_length_2d = 0.0

    # returns [
    #   [ [sketch_versions], [dep_node_id_1, cand_id_1], [dep_node_id_2, cand_id_2], ..., ],
    #   ...]
    def get_contexts_3d(self, stroke_id):
        #print("get_contexts_3d")
        # get different combinations from self.dependency_nodes[stroke_id].in_strokes
        # which belong to the same sketch_versions
        dep_node_id = self.stroke_id_to_dep_node_id[stroke_id]
        dep_node = self.dependency_nodes[dep_node_id]
        predecessor_dep_nodes_ids = [self.stroke_id_to_dep_node_id[i]
                                     for i in dep_node.in_strokes]
        predecessor_dep_nodes = [self.dependency_nodes[self.stroke_id_to_dep_node_id[i]]
                                 for i in dep_node.in_strokes]
        # initialize contexts_3d
        contexts_3d = [[cand_nb] for cand_nb in
                       predecessor_dep_nodes[0].candidate_nodes.keys()]
        #print("all_predecessors")
        #for i in range(len(predecessor_dep_nodes)):
        #    print("dep_node: "+str(i))
        #    curr_dep_node = predecessor_dep_nodes[i]
        #    print(curr_dep_node.is_assigned)
        #    print(curr_dep_node)
        #    for j in curr_dep_node.candidate_nodes.keys():
        #        curr_cand_node = curr_dep_node.candidate_nodes[j]
        #        print(curr_cand_node.sketch_versions)
        #print("contexts_3d: "+str(contexts_3d))
        for i in range(1, len(predecessor_dep_nodes)):
            new_contexts_3d = []
            curr_pred_node = predecessor_dep_nodes[i]
            # go through previous candidates and through all current candidates
            # and check for compatibility
            #print(i)
            #for cand_list in contexts_3d:
            #    print(cand_list)
            for cand_list in contexts_3d:
                found_any_candidates = False
                for curr_cand in curr_pred_node.candidate_nodes.values():
                    #print(curr_cand.sketch_versions)
                    # check with all candidates for cand_list
                    old_candidates = []
                    for cand_list_id, cand_id in enumerate(cand_list):
                        old_candidates.append(predecessor_dep_nodes[cand_list_id].
                                              candidate_nodes[cand_id])
                    old_candidates.append(curr_cand)
                    #print("old_candidates")
                    #for old_cand in old_candidates:
                    #    print(old_cand.sketch_versions)
                    # create set of sketch versions
                    compatible_versions = self.compute_compatible_versions(old_candidates)
                    if len(compatible_versions) == 0:
                        continue
                    #print("compatible_versions")
                    #print(compatible_versions)
                    found_any_candidates = True
                    new_cand_list = deepcopy(cand_list)
                    new_cand_list.append(curr_cand.candidate_nb)
                    new_contexts_3d.append(deepcopy(new_cand_list))
            contexts_3d = new_contexts_3d
        #for cand_list in contexts_3d:
        #    print(cand_list)

        #print("over")
        context_candidates = [[predecessor_dep_nodes[cand_id].candidate_nodes[cand]
                               for cand_id, cand in enumerate(context_3d)]
                              for context_3d in contexts_3d]
        context_versions = [self.compute_compatible_versions(context_candidate_list)
                            for context_candidate_list in context_candidates]
        #print(context_versions)
        contexts_3d = [Context3D(sketch_versions=context_versions[context_3d_id],
                                 dep_node_ids=predecessor_dep_nodes_ids,
                                 cand_node_ids=context_3d)
                       for context_3d_id, context_3d in enumerate(contexts_3d)]
        return contexts_3d

    def compute_compatible_versions(self, candidate_list):
        compatible_versions = set(candidate_list[0].sketch_versions)
        for cand in candidate_list[1:]:
            compatible_versions = set(
                cand.sketch_versions) & compatible_versions
        return list(compatible_versions)

    def display_all_candidates(self, dep_node_id):

        for cand_id in self.dependency_nodes[dep_node_id].candidate_nodes.keys():
            line = self.dependency_nodes[dep_node_id].candidate_nodes[cand_id].geometry
            sketch_3d = ps.register_curve_network("dep_node_id: "+str(dep_node_id)+", cand_id: "+str(cand_id),
                                                  nodes=np.array([line[0], line[-1]]),
                                                  edges=np.array([[0, 1]]),
                                                  enabled=False)

    def get_best_sketch(self):
        sketch_versions = self.get_n_best_sketch_versions(len(self.sketch_versions_reference.keys()))
        best_sketch_version = sketch_versions[0]
        lines = []
        intersections = []

        for dep_node_id, (cand_id, inter_set, score_container) in enumerate(self.sketch_versions_reference[best_sketch_version]):
            lines.append(self.dependency_nodes[dep_node_id].candidate_nodes[cand_id])
            intersections += inter_set
        return lines, intersections

    # display all sketch_versions
    def display_sketch_versions(self):
        # remove old displays
        ps.remove_all_structures()
        #for i in self.sketch_version_counter:
        # display in total_score order
        sketch_versions = self.get_n_best_sketch_versions(len(self.sketch_versions_reference.keys()))
        #sketch_versions = list(self.sketch_versions_reference.keys())[:10]
        #sketch_versions = self.get_n_best_sketch_versions(1)
        for sketch_version_id, sketch_version in enumerate(sketch_versions):
            # collect 3D lines
            lines = []
            scores = []
            line_coverages = []
            axis_alignments = []
            orthogonalities = []
            tangentialities = []
            planarities = []
            foreshortenings = []
            curve_geoms = []
            circularities = []
            is_assigned = []
            dep_node_ids = []

            for dep_node_id, (cand_id, _, score_container) in enumerate(self.sketch_versions_reference[sketch_version]):
                if len(self.dependency_nodes[dep_node_id].candidate_nodes) == 0:
                    continue
                dep_node_ids.append(dep_node_id)
                lines.append(self.dependency_nodes[dep_node_id].candidate_nodes[cand_id].geometry)
                #if dep_node_id == 26:
                #    print("dep_node_id: ", dep_node_id)
                #    print(lines[-1])
                scores.append(score_container.total_score)
                line_coverages.append(score_container.line_coverage)
                axis_alignments.append(score_container.axis_alignment)
                orthogonalities.append(score_container.orthogonality)
                tangentialities.append(score_container.tangentiality)
                planarities.append(score_container.planarity)
                foreshortenings.append(score_container.foreshortening)
                curve_geoms.append(score_container.curve_geom)
                circularities.append(score_container.circularity)
                is_assigned.append(self.dependency_nodes[dep_node_id].is_assigned)
            nodes = []
            edge_counter = 0
            edges = []
            enabled = False
            if sketch_version_id == 0:
                enabled = True

            line_ids = []
            for line_id, line in enumerate(lines):
                for p in line:
                    nodes.append(p)
                for p_id in range(len(line)-1):
                    edges.append([edge_counter, edge_counter+1])
                    line_ids.append(line_id)
                    edge_counter += 1
                edge_counter += 1

                #edge_counter += 1
                #nodes.append(line[0])
                #nodes.append(line[-1])
                #edges.append([edge_counter, edge_counter+1])
                #edge_counter += 2
            sketch_3d = ps.register_curve_network("sketch_version: "+str(sketch_version),
                                                  nodes=np.array(nodes), edges=np.array(edges),
                                                  enabled=enabled)
            sketch_3d.add_scalar_quantity("dep_node_ids", np.array(dep_node_ids)[line_ids], defined_on="edges",
                                          enabled=False, cmap="reds")
            sketch_3d.add_scalar_quantity("line_coverage", np.array(line_coverages)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("axis_alignment", np.array(axis_alignments)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("orthogonality", np.array(orthogonalities)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("tangentiality", np.array(tangentialities)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("planarity", np.array(planarities)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("foreshortening", np.array(foreshortenings)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("curve_geom", np.array(curve_geoms)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("circularity", np.array(circularities)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("is_assigned", np.array(is_assigned)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            sketch_3d.add_scalar_quantity("total_score", np.array(scores)[line_ids], defined_on="edges",
                                          enabled=True, cmap="jet", vminmax=(0., 1.))
            #print("sketch_version: "+str(sketch_version))
            #print("score: "+str(np.sum(scores)))
            #print(scores[-1])

            # collect 3D intersections
            points = []
            angles = []
            dep_node_ids = []
            distances = []
            curve_length = []
            for dep_node_id, (line_id, inter_set, _) in enumerate(self.sketch_versions_reference[sketch_version]):
                line = self.dependency_nodes[dep_node_id].candidate_nodes[line_id]
                line_length = tools_3d.line_3d_length(line.geometry)
                for inter in inter_set:
                    points.append(inter.coords_3d)
                    angles.append(inter.tangents_angle_3d)
                    dep_node_ids.append(self.stroke_id_to_dep_node_id[inter.stroke_ids])
                    dist = tools_3d.distance_point_to_polyline_vectorized(inter.coords_3d,
                                                               line.geometry)
                    distances.append(dist)
                    curve_length.append(line_length)
            if len(points) > 0:
                inter_cloud = ps.register_point_cloud("sketch_version: "+str(sketch_version),
                                        points=np.array(points), enabled=enabled, radius=0.01)
                inter_cloud.add_scalar_quantity("tangents_angle_3d", np.array(angles),
                                                enabled=True, cmap="jet", vminmax=(0., 90.))
                inter_cloud.add_scalar_quantity("first_stroke", np.array(dep_node_ids)[:, 0])
                inter_cloud.add_scalar_quantity("snd_stroke", np.array(dep_node_ids)[:, 1])
                inter_cloud.add_scalar_quantity("distance", np.array(distances))
                inter_cloud.add_scalar_quantity("line_length", np.array(curve_length))

    def __str__(self):
        output = "DependencyGraph"
        for dep_node_id, dep_node in enumerate(self.dependency_nodes):
            output += "\nDependencyNode: "+str(dep_node_id)
            output += str(dep_node)
        return output

    def update_dep_graph(self, merged_context_lines, dep_node_id):
        versions_used_by_stroke = []
        # go through all versions of all context_lines
        # if a version has already been used used by a previous context line,
        # create a new sketch_version

        # first: duplicate data when new branching occurs
        for context_line in merged_context_lines:
            for version_set_id, version_set in enumerate(context_line[2]):
                #print(version_set)
                for version in version_set:
                    if version in versions_used_by_stroke:
                        # create new version
                        self.sketch_version_counter += 1
                        self.sketch_versions_reference[self.sketch_version_counter] = \
                            deepcopy(self.sketch_versions_reference[version])
                        # update sketch_versions of candidate_nodes
                        for dep_node_id_tmp, (cand_id, _, _) in \
                                enumerate(self.sketch_versions_reference[version]):
                            self.dependency_nodes[dep_node_id_tmp].\
                                candidate_nodes[cand_id].sketch_versions.\
                                append(self.sketch_version_counter)
                        versions_used_by_stroke.append(self.sketch_version_counter)
                    else:
                        versions_used_by_stroke.append(version)

        #print("versions_used_by_stroke")
        #print(versions_used_by_stroke)
        # next: update data
        version_counter = 0
        for context_line_id, context_line in enumerate(merged_context_lines):
            versions_used_by_context_line = []
            for version_set_id, version_set in enumerate(context_line[2]):
                inter_set = context_line[1][version_set_id]
                for old_version in version_set:
                    version = versions_used_by_stroke[version_counter]
                    versions_used_by_context_line.append(version)
                    version_counter += 1
                    # update self.sketch_versions_reference
                    self.sketch_versions_reference[version].append([context_line_id, deepcopy(inter_set), ScoreContainer()])
                    # add intersections to prior strokes
                    for inter in inter_set:
                        prev_stroke_id = inter.stroke_ids[1 - np.argwhere(np.array(inter.stroke_ids) == self.dependency_nodes[dep_node_id].stroke_id).flatten()[0]]
                        self.sketch_versions_reference[version][self.stroke_id_to_dep_node_id[prev_stroke_id]][1].append(inter)
            cand = CandidateNode(candidate_nb=context_line_id)
            cand.geometry = context_line[0].geometry
            cand.plane_point = context_line[0].plane_point
            cand.plane_normal = context_line[0].plane_normal

            if (not self.dependency_nodes[dep_node_id].is_curve) and \
                    self.dependency_nodes[dep_node_id].axis_label < 3:
                # add axis-alignment score to candidate line
                cand.axis_alignment = tools_3d.compute_axis_alignment(cand.geometry,
                                                                      self.dependency_nodes[dep_node_id].axis_label)

            cand.sketch_versions = versions_used_by_context_line
            self.dependency_nodes[dep_node_id].insert_candidate_node(cand)
            self.stroke_lengths_3d[(dep_node_id, cand.candidate_nb)] = tools_3d.line_3d_length(cand.geometry)

        if len(list(self.dependency_nodes[dep_node_id].candidate_nodes.keys())) == 1:
            self.dependency_nodes[dep_node_id].is_assigned = True

    def remove_obsolete_sketch_versions(self):
        sketch_version_lengths = np.array([len(version) for version in self.sketch_versions_reference.values()])
        sketch_version_keys = np.array(list(self.sketch_versions_reference.keys()))
        curr_length = np.max(sketch_version_lengths)
        under_curr_length = np.argwhere(sketch_version_lengths < curr_length).flatten()
        #print("under_curr_length")
        #print(sketch_version_keys[under_curr_length])
        self.prune_sketch_versions(sketch_version_keys[under_curr_length])

    def update_stroke_score_function_curve(self, dep_node_id,
										   strokes_topology,
                                           extreme_intersection_distance,
                                           camera):

        max_nb_intersections = np.max([len(self.sketch_versions_reference[sketch_version][dep_node_id][1])
                                       for sketch_version in self.sketch_versions_reference.keys()])
        for sketch_version in self.sketch_versions_reference.keys():
            cand_id, inter_set, score_container = self.sketch_versions_reference[sketch_version][dep_node_id]
            cand_curve = self.dependency_nodes[dep_node_id].candidate_nodes[cand_id]

            if not self.dependency_nodes[dep_node_id].is_ellipse:
                # line coverage
                #inter_stroke_ids = [inter.stroke_ids for inter in inter_set]
                #print(inter_stroke_ids)
                stroke_id = self.dependency_nodes[dep_node_id].stroke_id
                arc_params = [inter.mid_inter_param[np.argwhere(np.array(inter.stroke_ids) == stroke_id).flatten()[0]]
                              for inter in inter_set]
                inter_dist = np.max(arc_params) - np.min(arc_params)
                line_coverage = 0.0
                if extreme_intersection_distance > 0.0:
                    # line_coverage can be greater than 1.0 if intersections from
                    # future strokes have been added
                    #print("sketch_version: ", sketch_version)
                    #print(dep_node_id)
                    #print(extreme_intersection_distance)
                    #print(arc_params)
                    #print(inter_dist)
                    line_coverage = min(1.0, inter_dist/extreme_intersection_distance)
                self.sketch_versions_reference[sketch_version][dep_node_id][2].line_coverage = line_coverage

            #nb_intersections = len(strokes_topology[
            #                           self.dependency_nodes[dep_node_id].stroke_id]["curve_intersections"])
            regular_intersections = [inter.tangents_angle_3d < 30.0 or inter.tangents_angle_3d > 80.0
                                     for inter in inter_set]
            geom_score = min(1.0, np.sum(regular_intersections)/max_nb_intersections)
            self.sketch_versions_reference[sketch_version][dep_node_id][2].curve_geom = geom_score

            if not self.dependency_nodes[dep_node_id].is_ellipse:
                view_dir = np.array(camera.view_dir)
                plane_normal = np.array(cand_curve.plane_normal)
                if np.isclose(np.linalg.norm(view_dir), 0.0) or \
                    np.isclose(np.linalg.norm(plane_normal), 0.0):
                    foreshortening_score = 0.0
                else:
                    view_dir /= np.linalg.norm(view_dir)
                    plane_normal /= np.linalg.norm(plane_normal)
                    foreshortening_score = np.abs(np.dot(view_dir, plane_normal))
                self.sketch_versions_reference[sketch_version][dep_node_id][2].foreshortening = foreshortening_score

            if self.dependency_nodes[dep_node_id].is_ellipse:
                plane_vecs = []
                plane_normal = cand_curve.plane_normal
                plane_normal /= np.linalg.norm(plane_normal)
                for i in range(3):
                    axis = np.zeros(3)
                    axis[i] = 1.0
                    if not np.isclose(1.0, np.dot(axis, plane_normal)):
                        plane_vec = np.cross(axis, plane_normal)
                        plane_vec /= np.linalg.norm(plane_vec)
                        plane_vecs.append(plane_vec)
                    if len(plane_vecs) == 2:
                        break

                circularity_score = 1.0 - tools_3d.get_ellipse_eccentricity(cand_curve.geometry,
                                                                            plane_vecs[0],
                                                                            plane_vecs[1])
                self.sketch_versions_reference[sketch_version][dep_node_id][2].circularity = circularity_score

            if self.dependency_nodes[dep_node_id].is_ellipse:
                total_score = 0.8*circularity_score + 0.2*geom_score
            else:
                total_score = 0.6*line_coverage + 0.2*geom_score + 0.2*foreshortening_score
            self.sketch_versions_reference[sketch_version][dep_node_id][-1].total_score = total_score

        if self.dependency_nodes[dep_node_id].is_ellipse:
            max_circularity = 0.0
            # normalize circularity score
            for sketch_version in self.sketch_versions_reference.keys():
                cand_id, inter_set, score_container = self.sketch_versions_reference[sketch_version][dep_node_id]
                max_circularity = max(max_circularity, score_container.circularity)
            for sketch_version in self.sketch_versions_reference.keys():
                self.sketch_versions_reference[sketch_version][dep_node_id][-1].circularity /= max_circularity
                self.sketch_versions_reference[sketch_version][dep_node_id][-1].total_score = \
                    0.8*self.sketch_versions_reference[sketch_version][dep_node_id][-1].circularity + \
                    0.2*self.sketch_versions_reference[sketch_version][dep_node_id][-1].curve_geom

    def update_stroke_score_function(self, dep_node_id, sketch, extreme_intersection_distance):
        for sketch_version in self.sketch_versions_reference.keys():
            cand_id, inter_set, score_container = self.sketch_versions_reference[sketch_version][dep_node_id]

            # line coverage
            arc_params = [inter.mid_inter_param[np.argwhere(inter.stroke_ids == self.dependency_nodes[dep_node_id].stroke_id).flatten()[0]]
                          for inter in inter_set]
            inter_dist = np.max(arc_params) - np.min(arc_params)
            line_coverage = 0.0
            if extreme_intersection_distance > 0.0:
                line_coverage = inter_dist/extreme_intersection_distance
            self.sketch_versions_reference[sketch_version][dep_node_id][2].line_coverage = line_coverage

            # axis_alignment
            axis_alignment = self.dependency_nodes[dep_node_id].candidate_nodes[cand_id].axis_alignment
            self.sketch_versions_reference[sketch_version][dep_node_id][2].axis_alignment = axis_alignment

            geom_score = axis_alignment
            if self.dependency_nodes[dep_node_id].axis_label == 3:
                # compute ortho_score
                ortho_score = 0.0
                tangentiality_score = 0.0
                planarity_score = 0.0

                for dep_inter in inter_set:
                    adj_inters = sketch.intersection_graph.get_adjacent_intersections(dep_inter.inter_id)
                    other_dep_node_ids = []
                    for adj_inter in adj_inters:
                        for i in adj_inter.stroke_ids:
                            if i != self.dependency_nodes[dep_node_id].stroke_id and \
                                    self.stroke_id_to_dep_node_id[i] < len(self.sketch_versions_reference[sketch_version]):
                                # check if this intersection has been used by other stroke
                                inter_ids_used = [inter.inter_id for inter in self.sketch_versions_reference[sketch_version][self.stroke_id_to_dep_node_id[i]][1]]
                                if adj_inter.inter_id in inter_ids_used:
                                    other_dep_node_ids.append(self.stroke_id_to_dep_node_id[i])
                    other_dep_node_ids = np.unique(other_dep_node_ids)
                    if len(other_dep_node_ids) == 0:
                        continue

                    other_cand_node_ids = [self.sketch_versions_reference[sketch_version][other_dep_node_id][0]
                                           for other_dep_node_id in other_dep_node_ids]

                    local_ortho_scores = [tools_3d.compute_gaussian_between_lines(
                        self.dependency_nodes[dep_node_id].candidate_nodes[cand_id].geometry,
                        self.dependency_nodes[other_dep_node_id].candidate_nodes[other_cand_id].geometry)
                        for (other_dep_node_id, other_cand_id) in zip(other_dep_node_ids, other_cand_node_ids)]
                    local_tangentiality_scores = [tools_3d.compute_gaussian_between_one_minus_lines(
                        self.dependency_nodes[dep_node_id].candidate_nodes[cand_id].geometry,
                        self.dependency_nodes[other_dep_node_id].candidate_nodes[other_cand_id].geometry)
                        for (other_dep_node_id, other_cand_id) in zip(other_dep_node_ids, other_cand_node_ids)]
                    if len(other_dep_node_ids) > 1:
                        local_planarity_scores = []
                        for comb in distinct_combinations(range(len(other_dep_node_ids)), 2):
                            plane_vec_1 = self.dependency_nodes[other_dep_node_ids[comb[0]]].candidate_nodes[other_cand_node_ids[comb[0]]].geometry
                            plane_vec_2 = self.dependency_nodes[other_dep_node_ids[comb[1]]].candidate_nodes[other_cand_node_ids[comb[1]]].geometry
                            local_planarity_scores.append(tools_3d.compute_planarity_score(
                                self.dependency_nodes[dep_node_id].candidate_nodes[cand_id].geometry,
                                plane_vec_1, plane_vec_2))
                        planarity_score = max(planarity_score, np.max(local_planarity_scores))
                    ortho_score = max(ortho_score, np.max(local_ortho_scores))
                    tangentiality_score = max(tangentiality_score, np.max(local_tangentiality_scores))

                geom_score = max(max(ortho_score, tangentiality_score), planarity_score)
                self.sketch_versions_reference[sketch_version][dep_node_id][2].orthogonality = ortho_score
                self.sketch_versions_reference[sketch_version][dep_node_id][2].tangentiality = tangentiality_score
                self.sketch_versions_reference[sketch_version][dep_node_id][2].planarity = planarity_score

            total_score = 0.4*line_coverage + 0.4*geom_score + 0.2*line_coverage*geom_score

            self.sketch_versions_reference[sketch_version][dep_node_id][-1].total_score = total_score

    # update scores for all strokes in dep_node_ids, i.e., only the strokes whose
    # score has been affected by the most recently added score
    def update_strokes_score_function(self, dep_node_ids, strokes_topology,
                                      extreme_intersections_distances_per_stroke,
                                      camera):
        for dep_node_id in dep_node_ids:
            if not self.dependency_nodes[dep_node_id].is_assigned:
                if not self.dependency_nodes[dep_node_id].is_curve:
                    self.update_stroke_score_function(dep_node_id, strokes_topology,
                                                      extreme_intersections_distances_per_stroke[dep_node_id])
                elif self.dependency_nodes[dep_node_id].is_curve or \
                        self.dependency_nodes[dep_node_id].is_ellipse:
                    self.update_stroke_score_function_curve(dep_node_id,
                                                            strokes_topology,
                                                            extreme_intersections_distances_per_stroke[dep_node_id],
                                                            camera)

        # update global score in the end
        for sketch_version in self.sketch_versions_reference.keys():
            scores = [self.dependency_nodes[dep_node_id].length_2d/self.total_length_2d * score_container.total_score
                      for dep_node_id, ( _, _, score_container) in
                      enumerate(self.sketch_versions_reference[sketch_version])]
            self.sketch_versions_scores[sketch_version] = np.sum(scores)

    def get_n_best_sketch_versions(self, n):
        keys = []
        scores = []
        for sketch_version in self.sketch_versions_reference.keys():
            scores.append(self.sketch_versions_scores[sketch_version])
            keys.append(sketch_version)
        best_indices = np.argsort(scores)
        keys = list(reversed(np.array(keys)[best_indices].tolist()))[:n]
        return keys

    def prune_sketch_versions(self, delete_sketch_versions):

        for sketch_version in delete_sketch_versions:
            for dep_node_id, (cand_id, _, score) in enumerate(self.sketch_versions_reference[sketch_version]):
                self.dependency_nodes[dep_node_id].candidate_nodes[cand_id].sketch_versions.remove(sketch_version)
            self.sketch_versions_reference.pop(sketch_version, None)
            self.sketch_versions_scores.pop(sketch_version, None)

        # remove empty candidate_nodes
        for dep_node_id in range(len(self.dependency_nodes)):
            for cand_id in list(self.dependency_nodes[dep_node_id].candidate_nodes.keys()):
                if len(self.dependency_nodes[dep_node_id].candidate_nodes[cand_id].sketch_versions) == 0:
                    self.dependency_nodes[dep_node_id].candidate_nodes.pop(cand_id, None)
                    self.stroke_lengths_3d.pop((dep_node_id, cand_id), None)
            if len(self.dependency_nodes[dep_node_id].candidate_nodes.keys()) == 1:
                self.dependency_nodes[dep_node_id].is_assigned = True

    # early_assignment takes as an input all impacted strokes by the recently
    # added stroke (so all intersected strokes). Its goal is to find if there
    # exists a set of candidate strokes which is clearly better than the rest.
    # We cluster all sketch_versions which have the same candidate_ids for
    # dep_node_ids strokes and take the max score for each cluster.
    # We then compare those scores using the Acceptance criteria from the paper.
    def early_assignment(self, dep_node_ids):
        ambiguous_dep_node_ids = []
        # don't consider already assigned dep_nodes
        for i in dep_node_ids:
            if not self.dependency_nodes[i].is_assigned:
                ambiguous_dep_node_ids.append(i)
        if len(ambiguous_dep_node_ids) < 2:
            return []
        # collect candidate_ids and scores for each sketch_version
        all_candidate_ids = np.zeros([len(self.sketch_versions_reference.keys()),
                                      len(ambiguous_dep_node_ids)], dtype=np.int)
        all_scores = np.zeros([len(self.sketch_versions_reference.keys()),
                                      len(ambiguous_dep_node_ids)])
        for sketch_version_id, sketch_version in enumerate(self.sketch_versions_reference.keys()):
            for vec_id, dep_node_id in enumerate(ambiguous_dep_node_ids):
                cand_id, _, score_container = self.sketch_versions_reference[sketch_version][dep_node_id]
                all_candidate_ids[sketch_version_id][vec_id] = cand_id
                all_scores[sketch_version_id][vec_id] = score_container.total_score
        # cluster sketch_versions based on candidate_ids
        all_sketch_versions = np.array(list(self.sketch_versions_reference.keys()))
        unique_candidate_ids, unique_candidate_ids_indices = np.unique(all_candidate_ids,
                                                                       axis=0,
                                                                       return_index=True)

        clustered_sketch_versions = []
        clustered_sketch_scores = []
        for unique_cand_id in unique_candidate_ids:
            cluster_ids = [vec_id for vec_id, cand_ids in enumerate(all_candidate_ids)
                           if np.all(np.equal(cand_ids, unique_cand_id))]
            clustered_sketch_versions.append(all_sketch_versions[cluster_ids])
            clustered_sketch_scores.append(np.max([np.mean(all_scores[cluster_id])
                                                  for cluster_id in cluster_ids]))

        sorted_scores = np.sort(clustered_sketch_scores)
        sorted_scores_ids = np.argsort(clustered_sketch_scores)
        accepted = acceptance_criteria(sorted_scores[-1], sorted_scores[-2])
        if accepted:
            # return all sketch_versions which do not belong to the best
            # sketch_version cluster: those sketch_versions will be pruned
            best_sketch_versions = clustered_sketch_versions[sorted_scores_ids[-1]]
            to_be_pruned_sketch_versions = [i for i in all_sketch_versions
                                            if not (i in best_sketch_versions)]
            return to_be_pruned_sketch_versions
        return []

    # accept_best_cand_line applies the acceptance_criteria to the most recently
    # added stroke.
    def accept_best_cand_line(self, dep_node_id):
        dep_node = self.dependency_nodes[dep_node_id]
        if dep_node.is_assigned:
            return []
        all_scores = []
        all_cand_ids = list(dep_node.candidate_nodes.keys())
        all_sketch_versions = list(self.sketch_versions_reference.keys())
        for cand_id in all_cand_ids:
            sketch_versions = dep_node.candidate_nodes[cand_id].sketch_versions
            max_score = np.max([self.sketch_versions_reference[version][dep_node_id][-1].total_score
                                for version in sketch_versions])
            all_scores.append(max_score)
        sorted_scores = np.sort(all_scores)
        sorted_scores_ids = np.argsort(all_scores)
        accepted = acceptance_criteria(sorted_scores[-1], sorted_scores[-2])
        if accepted:
            best_sketch_versions = dep_node.candidate_nodes[all_cand_ids[sorted_scores_ids[-1]]].sketch_versions
            to_be_pruned_sketch_versions = [i for i in self.sketch_versions_reference.keys()
                                            if i not in best_sketch_versions]
            return to_be_pruned_sketch_versions
        return []

def acceptance_criteria(q_best, q_snd_best):

    tau_score = 0.75
    tau_score_high = 0.98
    tau_ambiguity = 0.8

    #return (q_best > tau_score_high) or \
    #       ((q_best > tau_score) and (q_snd_best/q_best < tau_ambiguity))
    return ((q_best > tau_score_high) and (q_snd_best < tau_score_high)) or \
           ((q_best > tau_score) and (q_snd_best/q_best < tau_ambiguity))

def compute_normalized_dependency_graph(adj_mat, sketch):
    # find first intersection between two strokes with different axes
    # put those two strokes in the beginning
    # reshuffle adj_mat s.t. no stroke has no preceding intersection

    stroke_ordering = np.arange(adj_mat.shape[0])
    inter_stroke_pair = []
    for s_id in range(adj_mat.shape[0]):
        found_inter_stroke_pair = False
        intersecting_strokes = np.argwhere(adj_mat[s_id, :s_id]).flatten()
        for inter_s_id in intersecting_strokes:
            if sketch.strokes[s_id].axis_label != sketch.strokes[inter_s_id].axis_label and \
                sketch.strokes[s_id].axis_label < 3 and sketch.strokes[inter_s_id].axis_label < 3:
                found_inter_stroke_pair = True
                inter_stroke_pair = [inter_s_id, s_id]
                break
        if found_inter_stroke_pair:
            break
    adj_mat[[0, 1, inter_stroke_pair[0], inter_stroke_pair[1]]] = \
        adj_mat[[inter_stroke_pair[0], inter_stroke_pair[1], 0, 1]]
    stroke_ordering[[0, 1, inter_stroke_pair[0], inter_stroke_pair[1]]] = \
        stroke_ordering[[inter_stroke_pair[0], inter_stroke_pair[1], 0, 1]]

    # create reshuffled DependencyGraph
    dep_graph = DependencyGraph()
    # help_up_nodes are nodes which must be inserted after the last in_strokes
    # node was added to dep_graph
    help_up_nodes = []
    added_stroke_ids = []
    for vec_id, s_id in enumerate(stroke_ordering):
        neighbors = np.argwhere(adj_mat[s_id, :]).flatten()
        before_edges = neighbors[neighbors < s_id]
        after_edges = neighbors[neighbors > s_id]
        if vec_id > 0 and len(before_edges) == 0:
            if len(after_edges) == 0:
                # just ignore this stroke
                continue
            # get first intersection after the stroke has been drawn
            before_edges = [after_edges[0]]
            after_edges = after_edges[1:]
            help_up_nodes.append(DependencyNode(stroke_id=s_id,
                                                in_edges=before_edges,
                                                out_edges=after_edges,
                                                axis_label=sketch.strokes[s_id].axis_label))
        else:
            dep_graph.dependency_nodes.append(DependencyNode(stroke_id=s_id,
                                                             in_edges=before_edges,
                                                             out_edges=after_edges,
                                                             axis_label=sketch.strokes[s_id].axis_label))
            added_stroke_ids.append(s_id)
            # check if any held_up_node should be inserted now
            for h_node in reversed(help_up_nodes):
                strokes_added = [i in added_stroke_ids for i in h_node.in_strokes]
                if np.all(strokes_added):
                #if h_node.in_strokes[0] == s_id:
                    dep_graph.dependency_nodes.append(h_node)
                    added_stroke_ids.append(h_node.stroke_id)
                    #print("append: " + str(h_node.stroke_id))
                    help_up_nodes.remove(h_node)

    dep_graph.stroke_id_to_dep_node_id = np.zeros(len(sketch.strokes), dtype=np.int)
    for dep_node_id in range(len(dep_graph.dependency_nodes)):
        dep_graph.stroke_id_to_dep_node_id[dep_graph.dependency_nodes[dep_node_id].stroke_id] = dep_node_id
    return dep_graph

def compute_normalized_dependency_graph_v2(inter_graph, sketch, stroke_ordering):

    stroke_ordering = np.array(stroke_ordering)
    dep_graph = DependencyGraph()

    for vec_id, s_id in enumerate(stroke_ordering):
        neighbors = np.array(list(inter_graph.neighbors(s_id))).tolist()
        neighbors_ids = []
        for n in neighbors:
            n_id = list(np.argwhere(stroke_ordering == n).flatten())
            if len(n_id) > 0:
                neighbors_ids.append(n_id[0])

        before_edges = []
        after_edges = []
        for n_id in neighbors_ids:
            if n_id < vec_id:
                before_edges.append(stroke_ordering[n_id])
            else:
                after_edges.append(stroke_ordering[n_id])
        before_edges = np.array(before_edges)
        after_edges = np.array(after_edges)
        dep_graph.dependency_nodes.append(DependencyNode(stroke_id=s_id,
                                                         in_edges=before_edges,
                                                         out_edges=after_edges,
                                                         axis_label=sketch.strokes[s_id].axis_label))

    dep_graph.stroke_id_to_dep_node_id = np.zeros(len(sketch.strokes), dtype=np.int)
    for dep_node_id in range(len(dep_graph.dependency_nodes)):
        dep_graph.stroke_id_to_dep_node_id[dep_graph.dependency_nodes[dep_node_id].stroke_id] = dep_node_id
    return dep_graph

def init_dep_graph(strokes_topology):
    dep_graph = DependencyGraph()
    dep_graph.sketch_versions_reference[0] = []
    dep_graph.sketch_versions_scores[0] = []
    for s_id, s in enumerate(strokes_topology):
        if s["primitive_type"] != 0 or not s["depth_assigned"]:
            continue
        is_curve = False
        before_edges = []
        after_edges = []
        dep_node = DependencyNode(stroke_id=s_id,
                                  in_edges=before_edges,
                                  out_edges=after_edges,
                                  is_curve=is_curve)
        if s["depth_assigned"]:
            #print(s)
            dep_node.is_assigned = True
            candidate = CandidateNode(candidate_nb=0)
            candidate.geometry = s["stroke_3d"]
            candidate.sketch_versions = [0]
            dep_node.insert_candidate_node(candidate)
            dep_graph.sketch_versions_reference[0].append([0, [], ScoreContainer()])
            dep_graph.stroke_lengths_3d[(len(dep_graph.dependency_nodes), 0)] = tools_3d.line_3d_length(
                candidate.geometry)
        dep_graph.dependency_nodes.append(dep_node)

    scores = [score_container.total_score for (_, _, score_container) in
              dep_graph.sketch_versions_reference[0]]
    dep_graph.sketch_versions_scores[0] = np.sum(scores)
    dep_graph.stroke_id_to_dep_node_id = np.zeros(len(strokes_topology), dtype=np.int)
    for dep_node_id in range(len(dep_graph.dependency_nodes)):
        dep_graph.stroke_id_to_dep_node_id[dep_graph.dependency_nodes[dep_node_id].stroke_id] = dep_node_id

    return dep_graph

def append_to_dep_graph(dep_graph, stroke_ids, strokes_topology):

    for s_id in stroke_ids:
        s = strokes_topology[s_id]
        is_curve = False
        before_edges = []
        after_edges = []
        if s["primitive_type"] == 1:
            is_curve = True
            for inter in s["curve_intersections"]:
                before_edges.append(inter["straight_line_idx"])
            before_edges = np.unique(before_edges).tolist()
        dep_node = DependencyNode(stroke_id=s_id,
                                  in_edges=before_edges,
                                  out_edges=after_edges,
                                  is_curve=is_curve,
                                  is_ellipse=s["is_ellipse"])
        dep_graph.dependency_nodes.append(dep_node)
    for dep_node_id in range(len(dep_graph.dependency_nodes)):
        dep_graph.stroke_id_to_dep_node_id[dep_graph.dependency_nodes[dep_node_id].stroke_id] = dep_node_id
