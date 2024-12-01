function gs_data = readFields2D(filename)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Read output data in txt format and apply gather/scatter operations to
% have a continuous solution field.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear
clc

% Read the data into a matrix
data = readmatrix(filename);

% Extract columns into separate variables
x = data(:, 1);  % First column: x
y = data(:, 2);  % Second column: y
r = data(:, 3);  % Third column: r
u = data(:, 4);  % Fourth column: u
v = data(:, 5);  % Fifth column: v

% Combine x and y into a single matrix
xy = [x, y];

% Find unique (x, y) pairs and indices for grouping
[unique_xy, ~, group_indices] = unique(xy, 'rows');


% Compute averages for r, u, and v using accumarray
avg_r = accumarray(group_indices, r, [], @mean);
avg_u = accumarray(group_indices, u, [], @mean);
avg_v = accumarray(group_indices, v, [], @mean);

% Combine results into a single matrix for output
gs_data = [unique_xy, avg_r, avg_u, avg_v];


end

