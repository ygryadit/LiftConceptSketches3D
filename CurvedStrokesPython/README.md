# Install environment
```
conda env create -f curved_strokes_python.yml
```

# Input data

The input data has to be generated with the StraightLinesMatlab part of code.

# Launch command
```
python main.py -file_name=file_name.json
```

``` file_name.json ``` can be any reconstructed scaffold from ``` StraightLinesMatlab/folder_save```.

## Example
```
python main.py -file_name=../StraightStrokesMatLab/folder_save/student8/house/view1/student8_house_bestScore_full.json
```

# Output data
The output data will be saved in the ``` folder_save ``` folder.