% -------------------------
% CS5331: Aerial Computing
% Summer II, 2023
% -------------------------

function[] = proj_2(pp_strategy)


% network size: width (w) x height (h)
% e.g., 1000 (m) x 1000 (m) by default
w_begin= 0;
w_end = 1000;
h_begin = 0;
h_end = 1000;

strategy = pp_strategy;

% number of cells (subareas): 5-by-5, by default
n_cell = 5;
tot_cell = n_cell * n_cell;
size_cell = w_end / n_cell;

% n number of target points: 
n = 100;

% a set of rectangular subareas: 5-by-5
x_ = linspace(w_begin, w_end - size_cell, n_cell);
ux = [];
for i = 1:n_cell
    ux = [ux, x_]; 
end 
ux = ux';

y_ = ones(1, n_cell);  
uy = [];
for i = 1:n_cell
    uy = [uy, y_ .* (size_cell * (i - 1))];
end 
uy = uy';

% n number of weights: w, n-by-1, uniform
uw = ones(n, 1);

% n number of weights: w, n-by-1, uniform
% -- between the interval (w_begin, w_end) 
w_begin = 0;
w_end = 10;
w = w_begin + (w_end - w_begin) .* rand(n, 1);

% coverage area with radius, r (m), by default 100
r = 100;

% -----------------------
% scale-free distribution 
% -----------------------

% clustering exponent, alpha
alpha = 1.4;

% population, pop  
% -- initialize to zero
pop = ones(tot_cell, 1) - 1;

% probability, prob
% -- initialize to zero
prob = ones(tot_cell, 1) - 1;

% a set of rectangular subareas, 25-by-5
subarea_loc = [ux, uy];

% the first target point is randomly assigned to one of cells
pos_subarea = randi(tot_cell);

pos_x = randi(size_cell) + ux(pos_subarea);
pos_y = randi(size_cell) + uy(pos_subarea);
pop(pos_subarea) = pop(pos_subarea) + 1;

% the first target point - randomly assigned
loc(1, 1) = pos_x;
loc(1, 2) = pos_y;

% generate all scale-free target points (x, y)
for i = 2:n
    % calculate probabilities
    % -- sigma_pop = sum(pop, "all")
    sigma_pop = 0;
    for j = 1: tot_cell
        sigma_pop = sigma_pop + power(pop(j) + 1, alpha);
    end
    for j = 1: tot_cell
        prob(j) = power(pop(j) + 1, alpha) / sigma_pop; %power(sigma_pop, alpha);
        %prob(j) = power(pop(j), alpha) / power(sigma_pop, alpha)
    end
    % sanity check: if total probabilities are one
    %tot_prob = sum(prob, "all")

    % randomly choose one of subareas
    % -- pos_subarea = randi(tot_cell);
    
    % choose one of subareas based on the probability
    % -- generate a random and compare with cumulative probabilities 
    rand_prob = rand(1, 1); % generate between 0 to 1
    cumu_prob = 0; 
    for j = 1: tot_cell
        cumu_prob = cumu_prob + prob(j);
        if (cumu_prob >= rand_prob)
            pos_subarea = j;
            break
        end
    end

    % generate a position within the chosen subarea
    pos_x = randi(size_cell) + ux(pos_subarea);
    pos_y = randi(size_cell) + uy(pos_subarea);
    % increment the population of subarea
    pop(pos_subarea) = pop(pos_subarea) + 1;

    % add a new target point's (x, y) into a row
    loc = [loc; [pos_x, pos_y]];
end    

% loc2 is a matrix to store the counts of neighbors of each poi
loc2 = zeros(n,1);
for i = 1:n
    for j = 1+i : n
        d = sqrt((loc(i,1) - loc(j,1))^2 + (loc(i,2) - loc(j,2))^2);
        if d <= r && i~=j
            loc2(i) = loc2(i) + 1;
        end
    end
end

%loc3 matrix is used to map each poi with its number of neighbors
loc3 = zeros(size(loc2));
for i=1:n
    loc3(i,1)=i;
    loc3(i,2)=loc2(i);
end

%loc3 is sorted in descending order and stored in matrix loc2
loc4 = sortrows(loc3,2,'descend')
%matrix which has the final scan points in decreasing order of their
%neighbors
scanpoint = [0,0];

for i = 1:n
    temp = loc4(i,1);
    if temp ~= 0 %to skip the points whose value is set to 0
        figure(5);
        plot(loc(temp,1), loc(temp,2), "b^");
        hold on;
        scanpoint = [scanpoint; loc(temp,1), loc(temp,2)]; % adding to scanpoint matrix
                %plotting the circle around scan point
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + loc(temp,1);
        yunit = r * sin(th) + loc(temp,2);
        h = plot(xunit, yunit, "b");
        hold on;
        for j = 1 : n
            temp2 = loc4(j,1);
            if temp2 ~= 0
                d2 = sqrt((loc(temp,1) - loc(temp2,1))^2 + (loc(temp,2) - loc(temp2,2))^2);
                if d2 <= r && temp~=j
                    %the points inside the range of scan point are set to 0
                    loc4(j,:)=0;
                end   
            end
        end
    end
    loc2(i,:) = 0;
end
% draw target points
plot(loc(:, 1), loc(:, 2), "rx");
hold on;
set(gca, 'FontSize', 16);
set(gca, 'xTick', [-200:200:1200]);
set(gca, 'yTick', [-200:200:1200]);
xlabel('X', 'FontSize', 16);
ylabel('Y', 'FontSize', 16); 
axis equal
axis([-200 1200 -200 1200]);
title('ECP');
hold off;


num_points = length(scanpoint);

if strategy == 0
% Randomly select a point from the array to start the path. exclude the first point (0, 0) in the random selection as it is the initial point
    figure(6);
    plot(scanpoint(1,1), scanpoint(1,2), 'k-p', 'MarkerSize', 3, 'LineWidth', 1.5);
    hold on;
    random_indices = randperm(num_points - 1);

    path_x = [0];
    path_y = [0];
    current_x = 0;
    current_y = 0;
    total_distance = 0;
    for i = 1:num_points - 1
        next_point_index = random_indices(i);
        next_x = scanpoint(next_point_index,1);
        next_y = scanpoint(next_point_index,2);
        distance = sqrt((next_x - current_x)^2 + (next_y - current_y)^2);
        total_distance = total_distance + distance;
        plot(next_x, next_y, 'b^');
        hold on;
        path_x = [path_x, next_x];
        path_y = [path_y, next_y];
        line([path_x(i), next_x], [path_y(i), next_y], 'Color', 'r');
        hold on;
        current_x = next_x;
        current_x = next_y;
    end
fprintf('Total distance traveled by the drone in RAND: %.2f\n', total_distance);    
set(gca, 'FontSize', 16);
set(gca, 'xTick', [-200:200:1200]);
set(gca, 'yTick', [-200:200:1200]);
xlabel('X', 'FontSize', 16);
ylabel('Y', 'FontSize', 16); 
axis equal
axis([-200 1200 -200 1200]);
title('RAND');

hold off;


elseif strategy == 1
    total_distance = 0;
    figure(7);
    plot(scanpoint(1,1), scanpoint(1,2), 'k-p', 'MarkerSize', 3, 'LineWidth', 1.5);
    hold on;
    num_points = size(scanpoint, 1);
    path = [scanpoint(1,1), scanpoint(1,2)];
    % Creating a list to keep track of visited points
    visited = false(1, num_points);
    visited(1) = true; % Mark the initial point as visited
    % Find the nearest neighbor from the initial point (0, 0)
    current_point = 1;
    for i = 1:num_points-1
        min_distance = inf;
        nearest_neighbor = -1;    
        for j = 2:num_points
            if ~visited(j)
                distance = norm(scanpoint(current_point, :) - scanpoint(j, :));
                if distance < min_distance
                    min_distance = distance;
                    nearest_neighbor = j;
                end
            end
        end
    
    % Add the nearest neighbor to the path
        path = [path; scanpoint(nearest_neighbor, :)];
    
    % Mark the nearest neighbor as visited
        visited(nearest_neighbor) = true;
        total_distance = total_distance + min_distance;
    % Move to the nearest neighbor and continue finding the next nearest neighbor
        current_point = nearest_neighbor;
    end
plot(path(:, 1), path(:, 2), 'r-', 'LineWidth', 1);
hold on
plot(scanpoint(2:end, 1), scanpoint(2:end, 2), 'b^');

fprintf('Total distance traveled by the drone in NNF: %.2f\n', total_distance);
set(gca, 'FontSize', 16);
set(gca, 'xTick', [-200:200:1200]);
set(gca, 'yTick', [-200:200:1200]);
xlabel('X', 'FontSize', 16);
ylabel('Y', 'FontSize', 16); 
axis equal
axis([-200 1200 -200 1200]);
title('NNF');
hold off;
elseif strategy == 2
    total_distance = 0;
    num_points = size(scanpoint, 1);

% Initialize the path with the initial point (0, 0)
    path = [0, 0];
    figure(8);
    plot(scanpoint(1,1), scanpoint(1,2), 'k-p', 'MarkerSize', 3, 'LineWidth', 1.5);
    hold on;
    for i = 1:num_points
    % Get the coordinates of the current point
        current_x = scanpoint(i, 1);
        current_y = scanpoint(i, 2);
    
    % Plot the current point with a marker 'x' and color it red
        plot(current_x, current_y, 'b^');
        hold on;
    % Add the current point to the path
        path = [path; current_x, current_y];
    
    % If it's not the last point, draw a line connecting the current point to the next point
        if i < num_points
            next_x = scanpoint(i + 1, 1);
            next_y = scanpoint(i + 1, 2);
            line([current_x, next_x], [current_y, next_y], 'Color', 'r', 'LineWidth', 1);
            hold on;
        end
        distance = sqrt((next_x - current_x)^2 + (next_y - current_y)^2);
        total_distance = total_distance + distance;
    end
fprintf('Total distance traveled by the drone in DF: %.2f\n', total_distance);
set(gca, 'FontSize', 16);
set(gca, 'xTick', [-200:200:1200]);
set(gca, 'yTick', [-200:200:1200]);
xlabel('X', 'FontSize', 16);
ylabel('Y', 'FontSize', 16); 
axis equal
axis([-200 1200 -200 1200]);
title('DF');

hold off;
end



