% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% This code aim at analysing the variance of the GPMs
%% Include the GPML tool box and utils

cd gpml-matlab-v4.0-2016-10-19/
startup;
cd ..

addpath(genpath('utils'));



%% Simulate a random trajectory without noise

simulation_opt.imu_frequency = 100;
simulation_opt.duration = 2;
simulation_opt.data_overlap = 0.25;
simulation_opt.traj_profile = 'fast';
simulation_opt.one_axis = false;
simulation_opt.gravity_magnitude = 9.8;
simulation_opt.x_vel_offset = 2;
simulation_opt.acc_std = 0;
simulation_opt.gyr_std = 0;
simulation_opt.visualisation = false;


[ground_truth, imu_data] = SimulateTraj( simulation_opt );



%% Run the GPMs with different noise samples

acc_std = 0.02;
gyr_std = 0.002;

nb_run = 50;
rot_vel_pos = zeros(nb_run,9); 


quantum = 0.001;

for i = 1:nb_run
    
    fprintf(['.' num2str(i)]);
    
    imu_data_noisy = imu_data;
    imu_data_noisy.acc = imu_data_noisy.acc...
        + acc_std*randn(size(imu_data_noisy.acc));
    imu_data_noisy.gyr = imu_data_noisy.gyr...
        + gyr_std*randn(size(imu_data_noisy.gyr));
    
    
    gpm = Gpm(imu_data_noisy.acc,...
        imu_data_noisy.gyr,...
        imu_data_noisy.time,...
        ground_truth.start_time,...
        ground_truth.end_time,...
        quantum,...
        acc_std,...
        gyr_std,...
        simulation_opt.one_axis);
    
    
    rot_vel_pos(i,1:3) = LogMap(gpm.d_R);
    rot_vel_pos(i,4:6) = gpm.d_v';
    rot_vel_pos(i,7:9) = gpm.d_p';
end
disp(' ')
disp('Variances of GPM values around groundtuth')
disp(mat2str(var(rot_vel_pos,0,1)))
disp('GPM covariance diagonal')
disp(mat2str(diag(gpm.cov)))
