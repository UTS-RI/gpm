% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% Log mapping from SO(3) to angle-axis (vect 3x1)
function [ angle_axis ] = LogMap( rot_mat )

    trace_mat = trace(rot_mat);
    angle_axis = zeros(3,1);
    if trace_mat ~= 3
        phi = acos( (trace_mat - 1) / 2);
        skew_mat = (phi/(2*sin(phi))) * (rot_mat-rot_mat');
        angle_axis = [skew_mat(3,2,:);skew_mat(1,3,:);skew_mat(2,1,:)];
    end

end

