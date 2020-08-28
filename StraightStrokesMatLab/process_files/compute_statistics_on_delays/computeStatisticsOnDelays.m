function [percentage_immediate_decision,max_delay,adsps_mean,delays,hist_delays] = computeStatisticsOnDelays()

    folder_data = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v_48\v47planes_400';
    folder_data = fullfile(folder_data, 'sketches');

    designers = dir(folder_data);
    designers = {designers.name};
    designers = designers(3:end);

    
    max_delay = 0;
    num_same_stroke = 0;
    num_strokes = 0;
    adsps = [];
    hist_delays = [];
    delays = [];
    
    Num_assigned_in_the_end = [];
    num_strokes_all_skecthes = [];
    
    
    for d = 1:length(designers)
        objects = dir(fullfile(folder_data, designers{d}));
        objects = {objects.name};
        objects = objects(3:end);
        designers{d}
        for o = 1:length(objects )
            objects{o}
            folder_o = fullfile(folder_data, designers{d}, objects{o}, 'view1');
            filename = [designers{d} '_' objects{o}];
            
            if strcmp(designers{d}, 'Professional6') & (strcmp(objects{o}, 'tubes') | strcmp(objects{o}, 'wobble_surface'))
                continue;
            end
            
            try 
                [ strokes_topology, ~, ~] = readReconstructionJson( fullfile(folder_o, [filename '_confident_full.json'] ));
                [ strokes_topology_, ~, ~] = readReconstructionJson( fullfile(folder_o, [filename '_bestScore_full.json'] ));
            catch
                fullfile(folder_o, [filename '_bestScore_full.json'] )
                continue;
            end
            
            Num_assigned_in_the_end(end+1) = sum([strokes_topology_(:).depth_assigned]) - sum([strokes_topology(:).depth_assigned]);
            
            
            [max_delay_, num_same_stroke_, num_strokes_, adsps_, delays_, hist_delays] = computeDelaySketch(strokes_topology, hist_delays);
             if ~isempty(adsps_)
                adsps(end+1) = adsps_; 
             end
             if max_delay_(1) > max_delay(1)
                max_delay = max_delay_;
             end
             
             delays((end+1):(end+length(delays_))) = delays_;
             
             num_same_stroke = num_same_stroke + num_same_stroke_;
             num_strokes = num_strokes + num_strokes_;
%             num_strokes_all_skecthes(end+1) = num_strokes_ + Num_assigned_in_the_end(end);
             num_strokes_all_skecthes(end+1) = num_strokes_;
            
             percentage_immediate_decision = num_same_stroke/num_strokes;
%              max_delay
%              disp(adsps(end))
        end
    end
    
    percentage_immediate_decision = num_same_stroke/num_strokes;
    adsps_mean = mean(adsps);
end

function [max_delay, num_same_stroke, num_strokes, avg_delay_sketch_percentage, delays, hist_delays] = ...
                computeDelaySketch(strokes_topology, hist_delays)
            
    max_delay = -Inf;
    num_same_stroke = 0;
    num_strokes = 0;
    num_data_points = 0;
    delays = [];
    for i = 1:length(strokes_topology)
        if ~isempty(strokes_topology(i).created)
            num_strokes = num_strokes+1;
        end
        if ~isempty(strokes_topology(i).created) &  ~isempty(strokes_topology(i).assigned)
            
            i1 = strokes_topology(i).created;
            i2 = strokes_topology(i).assigned;
            delay = sum(([strokes_topology(i1:i2).primitive_type] == 0))-1
            strokes_topology(i).assigned - strokes_topology(i).created
%             delay = strokes_topology(i).assigned - strokes_topology(i).created;
           
            if delay > max_delay 
                max_delay(1) = delay;
                max_delay(2) = delay/length(strokes_topology);
               
            
            end
            
           
            if delay == 0
                num_same_stroke = num_same_stroke + 1;
            else
                num_data_points = num_data_points + 1;
                delays(num_data_points)  = delay;
                if delay > length(hist_delays)
                    hist_delays(delay) = 1;
                else
                    try
                        hist_delays(delay) = hist_delays(delay) + 1;
                    catch e
                        rethrow(e);
                    end
                end
             
            end
            
            
            

           
        end
    end
    
    avg_delay_sketch_percentage = mean(delays)/length(strokes_topology);
    if isempty(delays)
        avg_delay_sketch_percentage =[];
    end    
end