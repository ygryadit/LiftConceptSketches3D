clear all;
folder = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\backproject_reconstruction_json';
global folder_save;
global sketch_height;
sketch_height = 512;
foler_save_base = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\backproject_reconstruction_json';

[ strokes_topology, intersections, cam_param] = ...
readReconstructionJsonFelix( fullfile(folder, 'Professional2_vacuum_cleaner_bestScore_full.json' ));
[strokes_topology, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology, intersections, cam_param, []);
folder_save = fullfile(foler_save_base, 'Professional2_vacuum_cleaner');

% objectLFEffectAndSaveSVGFrames(strokes_topology, cam_param, 'lf')
objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology, cam_param, 'scaleXYZ')
% rotate3DAndSaveSVGFrames(strokes_topology, cam_param, '360');

[ strokes_topology, intersections, cam_param] = ...
readReconstructionJsonFelix( fullfile(folder, 'Professional6_mouse_bestScore_full.json' ));
[strokes_topology, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology, intersections, cam_param, []);
folder_save = fullfile(foler_save_base, 'Professional6_mouse');
% objectLFEffectAndSaveSVGFrames(strokes_topology, cam_param, 'lf')
objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology, cam_param, 'scaleXYZ')
% rotate3DAndSaveSVGFrames(strokes_topology, cam_param, '360');

[ strokes_topology, intersections, cam_param] = ...
readReconstructionJsonFelix( fullfile(folder, 'student8_house_bestScore_full.json' ));
[strokes_topology, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology, intersections, cam_param, []);
folder_save = fullfile(foler_save_base, 'student8_house');
%objectLFEffectAndSaveSVGFrames(strokes_topology, cam_param, 'lf')
objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology, cam_param, 'scaleXYZ')
% rotate3DAndSaveSVGFrames(strokes_topology, cam_param, '360');