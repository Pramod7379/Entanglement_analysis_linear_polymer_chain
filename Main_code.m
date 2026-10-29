% Script to run multiple scripts sequentially

% Define paths to your individual scripts

%  Execute the first script

scriptPaths = {
     'Master_code_first.m',
     };

    fprintf('Running script %d: %s\n', 1, scriptPaths{numel(scriptPaths)});
     run(scriptPaths{numel(scriptPaths)});  % Use 'run' to execute the script
    fprintf('Script %d completed. Created \output_mat_file\result_linku0.mat \n\n', 1);

%  Execute the second script

scriptPaths = {
     'Master_code_second.m',
     };

    fprintf('Running script %d: %s\n', 2, scriptPaths{numel(scriptPaths)});
     run(scriptPaths{numel(scriptPaths)});  % Use 'run' to execute the script
    fprintf('Script %d completed. Created \output_mat_file\seg_data.mat \n\n', 2);

%  Execute the third script

scriptPaths = {
     'Topological_entanglement.m',
     };

    fprintf('Running script %d: %s\n', 3, scriptPaths{numel(scriptPaths)});
     run(scriptPaths{numel(scriptPaths)});  % Use 'run' to execute the script
    fprintf('Script %d completed. Created \output_mat_file\topo.mat \n\n', 3);

%  Execute the fourth script

scriptPaths = {
     'Rheological_entanglement.m',
     };

    fprintf('Running script %d: %s\n', 4, scriptPaths{numel(scriptPaths)});
     run(scriptPaths{numel(scriptPaths)});  % Use 'run' to execute the script
    fprintf('Script %d completed. Created \output_mat_file\rheo.mat \n\n', 4);

fprintf('All scripts executed.\n');
