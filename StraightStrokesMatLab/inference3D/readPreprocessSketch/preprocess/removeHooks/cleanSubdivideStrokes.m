% function [strokes_cleaned] = cleanSubdivideStrokes(strokes, img)
% 
% This function removes hooks in the end of the strokes, it also might
% sibdivide the stroke if there are corners in the middle.

function [strokes_cleaned] = cleanSubdivideStrokes(strokes, img)
global SHOW_FIGS_PREPROCESS;
SHOW_FIGS_PREPROCESS_val  = SHOW_FIGS_PREPROCESS;
SHOW_FIGS_PREPROCESS = false;
if SHOW_FIGS_PREPROCESS
    hf1 = figure(3);
    hold off;
    imshow(img); 
    hold on;
else
    hf1 = [];
end

num_strokes = length(strokes);
j = 1;
imshow(img); 
hold on;
for i = 1:num_strokes 
    vertices = [cat(1,strokes(i).points(:).x), cat(1,strokes(i).points(:).y)];
    indices_org = 1:size(vertices,1);
    if isempty(vertices)
        continue;
    end
    [vertices, indices_sampled] = denslyResampleStroke(vertices);
    
    if size(vertices,1) < 3
        indices_set{1} = 1:size(vertices,1);
    else
%         fprintf('stroke %d\n', i)
        indices_set = cleanStroke(vertices, indices_org, indices_sampled, img, hf1);
    end
%     
%     figure(hf1);
%     plot(cat(1,strokes(i).points(:).x), cat(1,strokes(i).points(:).y), ...
%             'LineWidth', mean(cat(1,strokes(i).points(:).p)), 'Color', [0,0,1.0]);
%     
    
   
    for si = 1:length(indices_set)
        strokes_cleaned(j).points = strokes(i).points(indices_set{1});
        strokes_cleaned(j).mean_pressure = strokes(i).mean_pressure;
        strokes_cleaned(j).speed = strokes(i).speed;
        strokes_cleaned(j).primitive_type = strokes(i).primitive_type;
%         figure(hf);
%         plot(cat(1,strokes_cleaned(j).points(:).x), cat(1,strokes_cleaned(j).points(:).y), ...
%             'LineWidth', mean(cat(1,strokes_cleaned(j).points(:).p)), 'Color', [0,0,1.0]);
    
        j = j + 1;
    end
    
end

SHOW_FIGS_PREPROCESS = SHOW_FIGS_PREPROCESS_val;
end


function [vertices_new, ind_orignal] = denslyResampleStroke(vertices)
   vertices_new = [];
   ind_orignal =[];
    for i = 2:length(vertices)
        dir = vertices(i,:) - vertices(i-1,:);
        dist = norm(dir);
        k = ceil(dist/1.0);
        
        vertices_sampled = vertices(i-1,:) + [[0:k-1]'.*dir(1,1) [0:k-1]'.*dir(1,2)]/k;
        
        ind_orignal = [ind_orignal, size(vertices_new,1) + 1];
        vertices_new = [vertices_new; vertices_sampled];
%         figure(1);
%         hold off;
%         plot(vertices([i-1,i],1), vertices([i-1,2],2),'*-');
%         hold on;
%         plot(vertices_sampled(:,1), vertices_sampled(:,2),'o');
    end
    ind_orignal = [ind_orignal,  size(vertices_new,1) + 1];
    vertices_new(end+1,:) = vertices(end,:);  
end


function [indices_set] = cleanStroke(vertices, indices, indices_sampled, img, hf)
        
    global SHOW_FIGS_PREPROCESS;


    
%     k = polylineCurvature2D(vertices);
    if SHOW_FIGS_PREPROCESS
        figure(hf);
        hold off;
        imshow(img);
        hold on;
        plot(vertices(:,1), vertices(:,2), 'b');
    end
    [k, G0_points] = polylineCurvature2DFitCircleArrays(vertices, indices_sampled);
%     [k, G0_points] = polylineCurvature2DFitCircleArrays(vertices(indices_sampled,:), 1:length(indices_sampled));
    
    l_stroke = lengthStroke(vertices);

    if ~isempty(G0_points)
        max_k_ind = G0_points(1);
         if SHOW_FIGS_PREPROCESS
            plot(vertices(indices_sampled(G0_points),1), vertices(indices_sampled(G0_points),2), '*m');
         end
    else
         k = abs(k);
    
        curvature_sampled = zeros(size(vertices,1), 1);
        for jj = 1:length(indices_sampled)-1
            curvature_sampled(indices_sampled(jj):indices_sampled(jj+1)) = 0.5*(k(jj) + k(jj+1));
        end


        ind_cut = find( ( (k > 10*median(curvature_sampled)) & ...
                        (k > 0.125)));


        if SHOW_FIGS_PREPROCESS
            text(10, 10, sprintf('%.4f',median(curvature_sampled)));
            ind_cut_sampled = indices_sampled(ind_cut);

            for j = 1:length(ind_cut_sampled)
                plotCircle(vertices(ind_cut_sampled(j),1), vertices(ind_cut_sampled(j),2), 1.0/(k(ind_cut(j))), [0,0,0]);        
            end

            plot(vertices(ind_cut_sampled,1), vertices(ind_cut_sampled,2), '*');
            text(vertices(ind_cut_sampled,1), vertices(ind_cut_sampled,2), num2str(ind_cut));

           % Find max curvature:
          
        end
  
        %% Cut
        if isempty(ind_cut)
            % Check the angle in the beginning of the stroke:            
            vec1 = vertices(1,:) - vertices(2,:);
            vec2 = vertices(3,:) - vertices(2,:);
            vec1 = vec1/norm(vec1);
            vec2 = vec2/norm(vec2);
           
            if (dot(vec1, vec2) > 0)
                indices_set{1} = indices(2:end);
            else
                indices_set{1} = indices;
            end
            
            return;
        end

        cut_candidates = k(ind_cut);

        max_k_ind = ind_cut( ( cut_candidates == max(cut_candidates) ) );

        max_k_ind = max_k_ind(1);

        % Left and rigt strokes:
        if max_k_ind == 1
           max_k_ind = max_k_ind + 1;
        end

        if max_k_ind == length(indices_sampled)
           max_k_ind = max_k_ind - 1;
        end
    end
    
   
    
    
    
    left_indices_sampled = indices_sampled(1:max_k_ind);
    left_indices = indices(1:max_k_ind);
    right_indices_sampled = indices_sampled(max_k_ind:end) - indices_sampled(max_k_ind) + 1;
    right_indices = indices(max_k_ind:end);
      
    max_k_ind = indices_sampled(max_k_ind);
    
    
    left_vertices   = vertices(1:max_k_ind,:);
  
    right_vertices  = vertices(max_k_ind:end,:);
  
    
    % Check if a hook:
    l_left_stroke = lengthStroke(left_vertices);
    
    if (l_left_stroke > 0.15*l_stroke) 
        left_indices_set = cleanStroke(left_vertices, left_indices, left_indices_sampled, img, hf);
        
    else
        %hook
        left_indices_set = [];
    end
    
    l_right_stroke = lengthStroke(right_vertices);
    if (l_right_stroke > 0.15*l_stroke)
        right_indices_set = cleanStroke(right_vertices, right_indices, right_indices_sampled, img, hf);
      
    else
        %hook
        right_indices_set = [];
    end
    
    indices_set = [left_indices_set; right_indices_set];
end