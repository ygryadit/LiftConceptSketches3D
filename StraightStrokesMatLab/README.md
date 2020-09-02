# Input Data

The code is designed to run on the sketches from the [OpenSketch](https://repo-sam.inria.fr/d3/OpenSketch/) dataset.
You can download the data [here](https://repo-sam.inria.fr/d3/OpenSketch/Data/sketches/sketches_json_first_viewpoint.zip).
The zip archive has the follwoing structure:

```
sketches_json_first_viewpoint
└── designer_name_1
│   └── object_name_1
│   │   └── view1_concept.json
│    ...
│   └── object_name_m
│	└── view1_concept.json	
...
└── designer_name_k
``` 

We collected an additional set of sketches:  [OpenSketch++] (coming soon). The files in this archive are structured the same way as in the OpenSketch dataset.


# Launch command
To launch the code run in MatLab:

	launchOpenSketch(folder_sketches, designer_name, object_name, folder_save)


# Output Data
The code produces as the output the follwing folders and files:
```
folder_save
└── designer_name
    └── object_name
        └── view1
	    └── animation  			-- scg frames visualizing camera roration around the reconstructed shape
	    └── EstimatedVP
	    	└── lines_convergence.png 	-- An image that contains visualisation of three orthogonal directions
		└── vps.mat 	 		-- Mat file that contains the coordiantes of vanishing points.
	    └── images
		└── all_intersections_straight.svg     
		└── accuracy_radius.png 
		└── all_intersections_all_strokes.svg
		└── all_intersections_curved.svg  
		└── exact_intersections.png            
		└── intersections_likely.png           
		└── likely_intersections.svg           
		└── all_intersections_curved.svg       
	    └── lines_separation
	    	└── curves.svg  		-- this files is passed to a StrokeAggregator prior to curved strokes processing
		└── lines.svg   	    	
	    └── merged_svg
	    	└── lines_before_merging.svg  
	    	└── lines_after_merging.svg
	    └── svg_files
	    	└── lines_vp.svg		-- color coding of the lines for lines towards vanishign points and others
		
	    └── intersections_stat.mat  	-- contains fields: *num_intersections* and *num_likely_intersections*
	    └── preformance.mat  		-- *ellapsed time*
	    
	    └── designer_name_object_name_bestScore.obj                
	    └── designer_name_object_name_bestScore_full.json          
	    └── designer_name_object_name_bestScore_single_object.obj  
	    └── designer_name_object_name_confident.obj                
	    └── designer_name_object_name_confident_full.json          
	    └── designer_name_object_name_confident_single_object.obj  
	    └── designer_name_object_name_final.obj                    	-- the final reconstruction results used in the paper
	    └── designer_name_object_name_final_full.json               -- this json file is used as an input to the curve reconstruction part
	    └── designer_name_object_name_final_single_object.obj      
	    └── designer_name_object_name_highScore.obj                
	    └── designer_name_object_name_highScore_full.json          
	    └── designer_name_object_name_highScore_single_object.obj  	
```	    

## Files used as an input to the Curved strokes reconstruction part:

*folder_save/designer_name/object_name/view1/lines_separation/curves.svg*
and
*folder_save/designer_name/object_name/view1/designer_name_object_name_final_full.json*




# Code

## Folders structure

```
└── inference3D		-- the folder that contains all the functions responsible for lifting a sketch into 3D.
└── setup 		-- the folder that stores the scripts for setting up method constansts, parameters and flags.
```

## Visualization
To enable visualisation of each of the streps of the algorith change the following flags in *setDisplaySaveSettings.m* to *true*:
	
	global DISPLAY_INFO
	DISPLAY_INFO = false;

	global SHOW_FIGS
	SHOW_FIGS = false;

	global SHOW_FIGS_PREPROCESS
	SHOW_FIGS_PREPROCESS = false;


## Paper Section 4: Preprocessing
All the preprocessign described in Section 4 is done inside:

	[  	strokes_topology, ...
		intersections, ...
		cam_param ,...
		pairsInterInter,...
		ind_first_stroke] = ...
                	intialiseDataStructures();
			
			
where *strokes_topology*  and *intersections* are the two main data structures, described in detail below.
*cam_param* -- are the estimated camera parameters:
	
	cam_param = 

	  struct with fields:

			  P: [3×4 double] 	-- camera matrix
			  f: double		-- focal value	
	    principal_point: [3×1 double]
			  R: [3×3 double]	-- Rotation matrix
			  K: [3×3 double]	-- Calibration matrix
			  t: [3×1 double]	-- translation vector
		   view_dir: [3×1 double]
			fov: double		-- field of view in degrees
			  C: [3×1 double]	-- camera position
			 VP: [3×1 double]	-- vanishing points indices [1,2,3]
		   vp_coord: [3×2 double]	-- vanishing points coordinates
		   
*pairsInterInter* are the pairs of intersections that visually appear to be the same intersection, according to the criteria in the Supplemental Section 1.4.
*ind_first_stroke* an index of the stroke from which to start the reconstruction.

## Paper Section 5: Reconstructing Straight Strokes

### Core function    
 	[strokes_topology, intersections] = ...
                       assignDepthStroke(strokes_topology,....
                                         intersections,...
                                         ind_strk_zero_cnddts, ...
                                         cam_param,...
                                         pairsInterInter,...
                                         true);

# Main data structures:

## strokes_topology
	strokes_topology_ = 

	  num_strokes×1 struct array with fields:

	    accuracy_radius			: The radius for each stroke that defines the neighborhood in which the intersections visually 
	    					  group and the stroke close each other can be considered as intersecting
	    assigned				: the stoke number when the stroke is assigned
	    candidate_lines			: candidate_lines 
						  n_i×1 struct array with fields:
						    	origin: [x y z]
						    	dir: [x y z]
						    	coordinates3D_prior: [x1 y1 z1 x2 y2 z2]
						    	length3D: 
						    	configurations: [n×1 struct]
						    	inds_intrsctns:  All the intersections
							nds_intrsctns__assigned:  All the intersections with strokes with assigned 3D positions
							inds_intrsctns__mult_cnddts: All the intersections with strokes with candidate set of 3D positions*
							inds_intrsctns__mult_cnddts_ind: Version of the intersection (version of intersecting candidate line)
							inds_jnts_strks: List of strokes which contributed to computation of a joint cost
							p_intrsctns_dists: Distances from intersections to a proxy line
							p_directional: Reward associated with directional information
							p_coverage: Reward associated with 2D coverage of the strokes by the stroke piece between the first and last intersections
							p_full: Reward full of the strokes
							p_full_joint: Optimal joint reward
	
	    confidence				: measure of how much the best guess is better than the second best guess
	    created				: index of the stroke when the candaite lines for the current stroke were created for the first time
	    depth_assigned			: boolean flag *true* if stroke position is resolved
	    direction_vec			: the best directional vector that describes the stroke in 3D
	    ind_orth_ax				: Index of the axis [x, y, z] that is orthogonal to the plane in which the line lies
	    indcs_intrsctng_strks
	    indcs_intrsctns
	    inds_intrsctng_strks_eval
	    inds_intrsctns_eval
	    inter_seg_nums
	    length2D				: The length of polyline between the first and last intersections
	    length2DFull			: The length of a full 2D polyline
	    length2DLikely			: Porbably not used in the current version
	    length2DPrimitive			: Porbably not used in the current version	
	    length3D				: length in 3D
	    line_group				: Possible values:
						    	1: first vanishing point
							2: second vanishing point
							3: third vanishing point
							4: the rest of lines
							5: lines with a given Vanishing Point and a prior in which plane the line lies.
	    mean_pressure
	    merged_with
	    num_candidate_lines
	    num_candidate_lines_after_trim
	    num_candidate_lines_all
	    num_candidate_lines_before_trim
	    planes_normals
	    points2D				: num_pointsx1 struct array with fields:
							  p : pressure
							  t : time
							  x : x corrdinate
							  y : y coordinate
	    points2DOriginal
	    points3D
	    points3D_clean
	    poly2d_extended
	    primitive_geom			: The 2D coordinates of an approximating line geometry [x1 x2 y1 y2]
	    primitive_geom_3D
	    primitive_type			: Possible values:
							0: line
							1: curve
							2: mark
						       -1: grouped with some other stroke
						       -2: has zero intersections with any other stroke
						       -3: non defined speed

	    score
	    score_alignment
	    score_coverage
	    speed
    


## intersections

* .coordinates2D  
    *Intersection coordiantes*
* .strokes_indices
* .seg_nums
* .p_dist_str_segs
* .accuracy_radius
* .collinear
* .likely
* .coordinates3D
* .is_active
* .cnddts3D
    * .coordinates3D
        *Candidatee 3D coordinate*
    * .cnddt_lns 
        *2×1 cell array, first cell contains candidate lines indices from the first stroke in* .strokes_indices *and the second candidate lines from the second stroke*
    * .cnfgrtns 
        *2×1 cell array, each celll is again a cell* 


	
# Contact
If you have any questions please contact yulia.gryaditskaya@gmail.com
