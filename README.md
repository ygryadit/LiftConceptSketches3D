# LiftConceptSketches3D
The code for the paper: 
```
"Lifting Freehand Concept Sketches into 3D"
Yulia Gryaditskaya, Felix Hähnlein, Chenxi Liu, Alla Sheffer, Adrien Bousseau 
ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia), 2020.
```
# Description
The code consists of two parts the MatLab code for straight stroke reconstruction *StraightStrokesMatLab* and Python code for curved stroke reconstruction in *CurvedStrokesPython*.
See the README.md files for a detailed description of each of the parts.


# Project page
https://ns.inria.fr/d3/Lift3D/

# Input datasets
## OpenSketch

The code is designed to run on the sketches from the [OpenSketch](https://repo-sam.inria.fr/d3/OpenSketch/) dataset.
You can download the data [here](https://repo-sam.inria.fr/d3/OpenSketch/Data/sketches/sketches_json_first_viewpoint.zip).
The zip archive has the following structure:

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
## Additional design sketches OpenSketch++
We collected an additional set of sketches [OpenSketch++](https://repo-sam.inria.fr/d3/Lift3D/OpenSketch++.zip). The files in this archive are structured the same way as in the OpenSketch dataset.



# Reconstructed sketches are shown in the paper and the supplemental
All the reconstruction results are available [here](https://repo-sam.inria.fr/d3/Lift3D/reconstructions.zip).

# Newer version of the code
An update of the CurvedStrokesPython to support Python>=3.9 can be found in:
https://github.com/Liyb2002/LiftConceptSketches3D
