function [tr,envt]=get_trials(filename, num_trials)
    file_info = h5info(filename);
    group_name = file_info.Groups.Name;
    tr = h5read(filename, [group_name,'/trials']);
    tr=tr(1:2,1:num_trials);
 %to get the start and end states of the environment for each trial, we do:
    envt = h5read(filename, [group_name,'/trial_info'], [1 1], [2 Inf]);
    % row 1 of the cell corresponds to starting state of the environment
    % row 2 of the cell corresponds to ending state of the environment
    % a value of 1 indicates H+, -1 indicates H-
    
  % example of extracting click times for left and right streams (col vecs)
    % for trial 50
    % [left_train, right_train] = tr{:,50};
end