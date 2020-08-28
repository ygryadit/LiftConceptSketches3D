# LiftConceptSketches3D
The code for the paper: 
```
"Lifting Freehand Concept Sketches into 3D"
Yulia Gryaditskaya, Felix Hähnlein, Chenxi Liu, Alla Sheffer, Adrien Bousseau 
ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia), 2020.
```
# Description
The code consists of two parts the MatLab code for straight strokes reconstruction *StraightStrokesMatLab* and Python code for curved stroke recosntruction in *CurvedStrokesPython*.
See the README.md files for detailed description of each of the parts.

We are currently working on the unified implementation in Python, which we hope to release in the future.

# Input datasets
## OpenSketch

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
## Additional design sketches OpenSketch++
We additionally collected a set of sketches, which are available [coming soon]. The files in this archive are structured the same way as in the OpenSketch dataset.



# Reconstructed sketches shown in the paper and the supplemental

