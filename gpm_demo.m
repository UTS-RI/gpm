% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% This code aims at demonstrating the computation of GPMs
%% Include the GPML tool box and utils

cd gpml-matlab-v4.0-2016-10-19/
startup;
cd ..

addpath(genpath('utils'));



%% Simulate a random trajectory

simulation_opt.imu_frequency = 100;
simulation_opt.duration = 2;
simulation_opt.data_overlap = 0.25;
simulation_opt.traj_profile = 'fast';
simulation_opt.one_axis = false;
simulation_opt.gravity_magnitude = 9.8;
simulation_opt.x_vel_offset = 2;
simulation_opt.acc_std = 0.02;
simulation_opt.gyr_std = 0.002;
simulation_opt.visualisation = true;


[ground_truth, imu_data] = SimulateTraj( simulation_opt );


%% Compute the GPM and PM

pm = Pm(imu_data.acc,...
        imu_data.gyr,...
        imu_data.time,...
        ground_truth.start_time,...
        ground_truth.end_time,...
        simulation_opt.acc_std,...
        simulation_opt.gyr_std);


% Quantum of time for numerical integration of the rotation in the
% multi-axis rotation case   
quantum = 0.001;

gpm = Gpm(imu_data.acc,...
        imu_data.gyr,...
        imu_data.time,...
        ground_truth.start_time,...
        ground_truth.end_time,...
        quantum,...
        simulation_opt.acc_std,...
        simulation_opt.gyr_std,...
        simulation_opt.one_axis);

    
%% Compute and display errors
    
gpm_pos_error = norm(gpm.d_p - ground_truth.d_p);
gpm_rot_error = norm(LogMap(gpm.d_R'*ground_truth.d_R));

pm_pos_error = norm(pm.d_p - ground_truth.d_p);
pm_rot_error = norm(LogMap(pm.d_R'*ground_truth.d_R));

disp(['PM error: ' num2str(pm_pos_error) ' m    '...
    num2str(pm_rot_error*180/pi) ' deg']);
disp(['GPM error: ' num2str(gpm_pos_error) ' m    '...
    num2str(gpm_rot_error*180/pi) ' deg']);
