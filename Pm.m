% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% Function to compute the standard preintegration
function [ pm ] = Pm(  acc, gyr, time, start_time, query_time,...
    acc_sd, gyr_sd)


    % Mask the data to get only the bit between the start and query time
    mask_time = (time>=start_time) & (time<query_time);
    t = time( mask_time );
    t = [t; query_time];
    if t(1) == start_time
        mask_data = mask_time;
        d_t = t(2:end) - t(1:(end-1));
    else
        mask_data = [mask_time(2:end); false] | mask_time;
        d_t = t - [start_time; t(1:(end-1))];
    end
    acc_raw = acc(mask_data, :);
    gyr_raw = gyr(mask_data, :);
    acc_bias = zeros(3,1);
    gyr_bias = zeros(3,1);
    
    
     % Initialise the output
    d_R = eye(3);
    d_v = zeros(1,3);
    d_p = zeros(1,3);
    cov_pre_int = zeros(9);
    
    
    acc_cov = acc_sd^2 * ones(3,1);
    gyr_cov = gyr_sd^2 * ones(3,1);
    
    
    % Iterative preintegration
    for i = 1:length(d_t)
        acc_rot = d_R*((acc_raw(i,:))' - acc_bias);
        d_p = d_p...
            + (d_v * d_t(i))...
            + (acc_rot' * d_t(i) * d_t(i))/2;
        d_v = d_v + (acc_rot' * d_t(i));
        
        e_R = ExpMap((gyr_raw(i,:) - gyr_bias') * d_t(i));
        j_r = So3RightJacobian((gyr_raw(i,:) - gyr_bias') * d_t(i));
        d_R_ik = d_R;
        d_R = d_R * e_R;
        
        
        % Covariance part
        cov_imu = [diag(gyr_cov), zeros(3);...
            zeros(3), diag(acc_cov)];
        
        skew_acc = [0, -acc_raw(i,3), acc_raw(i,2);...
                    acc_raw(i,3), 0, -acc_raw(i,1);...
                   -acc_raw(i,2), acc_raw(i,1), 0];
        
        A = [e_R', zeros(3), zeros(3);...
            -d_R_ik*skew_acc*d_t(i), eye(3), zeros(3);...
            -d_R_ik*skew_acc*d_t(i)*d_t(i)/2 , eye(3), eye(3)];
        
        B = [j_r*d_t(i), zeros(3);...
            zeros(3), d_R_ik*d_t(i);...
            zeros(3), d_R_ik*d_t(i)*d_t(i)/2];
            
        cov_pre_int...
            = (A * cov_pre_int * A') + (B * cov_imu * B');
    end
    
    % Fill output
    pm.d_R = d_R;
    pm.d_v = d_v';
    pm.d_p = d_p';
    pm.cov = cov_pre_int;

end


function j_r = So3RightJacobian(r_vect)
    j_r = eye(3);
    r_norm = norm(r_vect);
    if r_norm ~= 0
        r_skew = [0, -r_vect(3), r_vect(2);...
                  r_vect(3), 0, -r_vect(1);...
                 -r_vect(2), r_vect(1), 0]; 
        j_r = j_r - ((1-cos(r_norm))/(r_norm^2))*r_skew...
            + ((r_norm - sin(r_norm))/(r_norm^3))*r_skew*r_skew;
    end
end
