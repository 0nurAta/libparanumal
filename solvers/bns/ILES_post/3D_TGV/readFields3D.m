function gs_data = readFields3D(filename)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Read output data in txt format and apply gather/scatter operations to
% have a continuous solution field.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Read the data into a matrix
data = readmatrix(filename);

% Extract columns into separate variables
x = data(:, 1);  % First column: x
y = data(:, 2);  % Second column: y
z = data(:, 3);  % Third column: z
r = data(:, 4);  % Fourth column: u
u = data(:, 5);  % Fifth column: v
v = data(:, 6);  
w = data(:, 7);
% Combine x, y and z into a single matrix
xyz = [x, y, z];

% Find unique (x, y) pairs and indices for grouping
[unique_xyz, ~, group_indices] = unique(xyz, 'rows');


% Compute averages for r, u, and v using accumarray
% avg_r = accumarray(group_indices, r, [], @mean);
avg_u = accumarray(group_indices, u, [], @mean);
avg_v = accumarray(group_indices, v, [], @mean);
avg_w = accumarray(group_indices, w, [], @mean);

% Combine results into a single matrix for output
gs_data = [unique_xyz, avg_u, avg_v, avg_w];

    % Write the output to a text file
    output_filename = "a.txt";
    fileID = fopen(output_filename, 'w');
    if fileID == -1
        error('Unable to open file %s for writing.', output_filename);
    end

    % Write a header
    fprintf(fileID, 'unique_xyz, avg_u, avg_v, avg_w\n');

    % Write the data row by row
    fprintf(fileID, '%.6f %.6f %.6f %.6f %.6f %.6f\n', gs_data');

    % Close the file
    fclose(fileID);

    fprintf('Processed data has been saved to %s\n', output_filename);
end

