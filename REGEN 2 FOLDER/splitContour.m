%% Contour Splitter
% Author: Andrew Radulovich (aradulov@purdue.edu), Kamon Blong
% First Created: 04/14/2023
% Last Updated: 01/21/2024

%{ 

Description: Takes in engine geometry with a variable distance between
    points and converts it to a contour with equal-length segments.

%}

function [x_contour_new, r_contour_new, L_seg] = splitContour(x_contour, r_contour, x_chamber, x_converging, x_diverging, y_chamber, y_converging,y_diverging, resolution, arcbool, debug)

%% Initialize Variables

num_points_chamber = resolution(1); % number of equally spaced points= 
num_points_converging = resolution(2);
num_points_diverging = resolution(3);

%% Initialize points

% remove repeated indices
[~, unique_index] = unique(x_converging);
unique_index_sorted = sortrows(unique_index);
x_converging = x_converging(unique_index_sorted);
y_converging = y_converging(unique_index_sorted);

% remove repeated indices
[~, unique_index] = unique(x_diverging);
unique_index_sorted = sortrows(unique_index);
x_diverging = x_diverging(unique_index_sorted);
y_diverging = y_diverging(unique_index_sorted);


%% Calculations
if arcbool 

% calculate the cumulative distance for each point
    cumulative_distance = zeros(1, length(x_contour));
    for i = 2:length(x_contour)
        cumulative_distance(i) = cumulative_distance(i-1) + sqrt((x_contour(i)-x_contour(i-1))^2 + (r_contour(i)-r_contour(i-1))^2);
    end
    
    % Create an interpolation function for x and r points based on cumulative distance
    x_interp = interp1(cumulative_distance, x_contour, 'linear', 'pp');
    
    r_interp = interp1(cumulative_distance, r_contour, 'linear', 'pp');
    
    % calculate the new cumulative distances for the equally spaced points
    L_seg = cumulative_distance(end) / (num_points - 1);
    new_cumulative_distance = 0:L_seg:cumulative_distance(end);
    
    % interpolate the new x and r points using the new cumulative distances
    x_contour_new = ppval(x_interp, new_cumulative_distance);
    r_contour_new = ppval(r_interp, new_cumulative_distance);
else 
    L_chamber = abs(x_chamber(end) - x_chamber(1));
    L_converging = abs(x_converging(end) - x_converging(1));
    L_diverging = abs(x_diverging(end) - x_diverging(1));
    L_seg(1) = (L_chamber/(num_points_chamber-1)); % change in distance per step [m]
    L_seg(2) = (L_converging/(num_points_converging-1)); % change in distance per step [m]
    L_seg(3) = (L_diverging/(num_points_diverging-1)); % change in distance per step [m]
    x_chamber_new = x_chamber(1):L_seg(1):x_chamber(end); % length vector
    x_converging_new = x_converging(1):L_seg(2):x_converging(end); % length vector
    x_diverging_new = x_diverging(1):L_seg(3):x_diverging(end); % length vector
    r_chamber_new = interp1(x_chamber,y_chamber,x_chamber_new,'linear'); % linearly interpolate radius vector 
    r_converging_new = interp1(x_converging,y_converging,x_converging_new,'linear'); % linearly interpolate radius vector 
    r_diverging_new = interp1(x_diverging,y_diverging,x_diverging_new,'linear'); % linearly interpolate radius vector 
    
    % combine all parts
    x_total_new = [x_chamber_new, x_converging_new, x_diverging_new];
    r_total_new = [r_chamber_new, r_converging_new, r_diverging_new];

    % remove repeated indices
    [~, unique_index] = unique(x_total_new);
    unique_index_sorted = sortrows(unique_index);
    x_contour_new = x_total_new(unique_index_sorted);
    r_contour_new = r_total_new(unique_index_sorted);

end

%% Display Outputs

if debug
    % print the equal spacing distance
    disp(['Equal spacing distance: ', num2str(L_seg)]);
    
    % plot the original and new points
    figure;
    plot(x_contour, r_contour, 'bo-', 'LineWidth', 1);
    hold on;
    plot(x_contour_new, r_contour_new, 'ro-', 'LineWidth', 1);
    legend('Original Points', 'Equally Spaced Points');
    xlabel('Inches X');
    ylabel('Inches Y');
    axis equal;
    grid on;
    title('Engine Contour with Equally Spaced Points');
end