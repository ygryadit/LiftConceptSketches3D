# Input Data

The code is designed to run on the sketches from the [OpenSkecth](https://repo-sam.inria.fr/d3/OpenSketch/) dataset.
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

We additionally collected a set of sketches, which are available [coming soon]. The files in this archive are structured the same way as in the OpenSketch dataset.


# Launch command
To launch the code run in MatLab:

	launchOpenSketch(folder_skecthes, designer_name, object_name, folder_save)


# Output Data
The code produces as the output the follwing folders and files:
```
folder_save
└── designer_name
    └── object_name
        └── view1
	    └── EstimatedVP
	    	└── lines_convergence.png 	-- An image that contains visualisation of three orthogonal directions
		└── vps.mat 	 		-- Mat file that contains the coordiantes of vanishing points.
	    └── images
	    └── lines_separation
	    └── merged_svg
	    └── svg_files
```	    





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
All the preprocessign decribed in Sectino 4 is done inside:

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
			  f: double		-- foval value	
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


## Core function    
    [strokes_topology, intersections, numbers_of_strokes_to_revisit] = ...
               assignDepthtoStrokePrecisely(strokes_topology,....
                                            intersections,...
                                            stroke_ind_cur, ...
                                            sketch_height, ...
                                            sketch_width,...
                                            cam_param,...
                                            img,...
                                            folder_save, ...
                                            pairsInterInter,...
                                            do_assign_depth)

# Main data structures:

## strokes_topology
    

* *points2D*
    * .p -- pressure
    * .t -- time
    * .x -- x corrdinate
    * .y -- y coordinate
    
* *primitive_type*      
    * 0: line
    * 1: curve
    * 2: mark
	* -1: grouped with some other stroke
	* -2: has zero intersections with any other stroke
	* -3: non defined speed
   
* *primitive_geom*
    * In case of a line, the 2D coordinates of an approximating line geometry [x1 x2 y1 y2]
        
* *mean_pressure*
* *length2DFull* 
    *   The length of a full 2D polyline 
*  *length2D* 
    *   The length of polyline between the first and last intersections
*  *length3D*
*  *line_group*
    * 1: first vanishing point
    * 2: second vanishing point
    * 3: third vanishing point
    * 4: the rest of lines
	* 5: lines with a given Vanishing Point and a prior in which plane the line lies.
* *speed*
    *  An average speed of a stroke
* *accuracy_radius*
    *  The radius for each stroke that defines the neighborhood in which the intersections visually group and the stroke close each other can be considered as intersecting
* indcs_intrsctng_strks
* indcs_intrsctns
* inds_intrsctns_eval
* inds_intrsctng_strks_eval
* primitive_geom_3D
* points3D
* score
* score_alignment
* score_coverage
* confidence
* direction_vec
* depth_assigned
* *.candidate_lines* [n×1 struct]
    * .origin: [x y z]
    * .dir: [x y z]
    * .coordinates3D_prior: [x1 y1 z1 x2 y2 z2]
    * .length3D: 
    * .configurations: [n×1 struct]
        *   .inds_intrsctns:                    
        *All the intersections*
        *   .inds_intrsctns__assigned:          // *All the intersections with strokes with assigned 3D positions*
        *   .inds_intrsctns__mult_cnddts:       // *All the intersections with strokes with candidate set of 3D positions*
        *   .inds_intrsctns__mult_cnddts_ind    // version of the intersection (version of intersecting candidate line)
        *   .inds_jnts_strks:                   // list of strokes which contributed to computation of a joint cost
        *   .p_intrsctns_dists:                 // distances from intersections to a proxy line
        *   .p_directional:                     // reward associated with directional information
        *   .p_coverage:                        // reward associated with 2D coverage of the strokes by the stroke piece between the first and last intersections
        *   .p_full:                            // reward full of the strokes
        *   .p_full_joint:                      // optimal joint reward
* *ind_orth_ax*
	* Index of the axis [x, y, z] that is orthogonal to the plane in which the line lies.
	
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


	
