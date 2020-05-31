% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% Exponential mapping from axis-angle (vect 3x1) to SO(3)
function [ rot_mat ] = ExpMap( angle_axis )
    
    rot_mat = eye(3);
    norm_vect = norm(angle_axis);
    if norm_vect~=0
        sMat = [0, (-angle_axis(3)), angle_axis(2);...
                angle_axis(3), 0, (-angle_axis(1));...
                (-angle_axis(2)), angle_axis(1), 0];
        rot_mat = rot_mat...
            + ( (sin(norm_vect)/norm_vect) * sMat)...
            + ( ( (1-cos(norm_vect)) / (norm_vect^2)) * sMat * sMat);
    end
        
    
end

