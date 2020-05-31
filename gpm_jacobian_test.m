% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% This code aims at analysing the postintegration correction Jacobians
% computed with the GPMs
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
simulation_opt.visualisation = false;


[ground_truth, imu_data] = SimulateTraj( simulation_opt );


%% Compute numerical jacobians

diff_delta = 0.001;


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

    
% Time-shift Jacobian
imu_data_shift = imu_data;
gpm_dt = Gpm(imu_data_shift.acc,...
        imu_data_shift.gyr,...
        imu_data_shift.time,...
        ground_truth.start_time + diff_delta,...
        ground_truth.end_time + diff_delta,...
        quantum,...
        simulation_opt.acc_std,...
        simulation_opt.gyr_std,...
        simulation_opt.one_axis);
delta_d_R_d_t_num = LogMap(gpm.d_R'*gpm_dt.d_R)/diff_delta;
delta_d_v_d_t_num = (gpm_dt.d_v - gpm.d_v)/diff_delta;
delta_d_p_d_t_num = (gpm_dt.d_p - gpm.d_p)/diff_delta;


disp('Jacobian d_R time-shift')
disp(['Analytical  ' mat2str(gpm.delta_d_R_d_t)])
disp(['Numerical   ' mat2str(delta_d_R_d_t_num)])
disp(' ')
disp('Jacobian d_v time-shift')
disp(['Analytical  ' mat2str(gpm.delta_d_v_d_t)])
disp(['Numerical   ' mat2str(delta_d_v_d_t_num)])
disp(' ')
disp('Jacobian d_p time-shift')
disp(['Analytical  ' mat2str(gpm.delta_d_p_d_t)])
disp(['Numerical   ' mat2str(delta_d_p_d_t_num)])

% Gyr bias Jacobian
gpm_dbw = cell(3,1);
delta_d_R_d_bw_num = zeros(3);
delta_d_v_d_bw_num = zeros(3);
delta_d_p_d_bw_num = zeros(3);
for i = 1:3
    if (~simulation_opt.one_axis)...
            || ( (simulation_opt.one_axis) && (i == 3) )
        imu_data_shift = imu_data;
        imu_data_shift.gyr(:,i) = imu_data_shift.gyr(:,i) + diff_delta;
        gpm_dbw{i} = Gpm(imu_data_shift.acc,...
                imu_data_shift.gyr,...
                imu_data_shift.time,...
                ground_truth.start_time,...
                ground_truth.end_time,...
                quantum,...
                simulation_opt.acc_std,...
                simulation_opt.gyr_std,...
                simulation_opt.one_axis);
        delta_d_R_d_bw_num(:,i) = LogMap(gpm.d_R'*gpm_dbw{i}.d_R)/diff_delta;
        delta_d_v_d_bw_num(:,i) = (gpm_dbw{i}.d_v - gpm.d_v)/diff_delta;
        delta_d_p_d_bw_num(:,i) = (gpm_dbw{i}.d_p - gpm.d_p)/diff_delta;
    end
end

disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('Jacobian d_R gyr bias')
disp('Analytical  ')
disp(mat2str(gpm.delta_d_R_d_bw(1,:)))
disp(mat2str(gpm.delta_d_R_d_bw(2,:)))
disp(mat2str(gpm.delta_d_R_d_bw(3,:)))
disp('Numerical   ')
disp(mat2str(delta_d_R_d_bw_num(1,:)))
disp(mat2str(delta_d_R_d_bw_num(2,:)))
disp(mat2str(delta_d_R_d_bw_num(3,:)))
disp(' ')
disp('Jacobian d_v gyr bias')
disp('Analytical  ')
disp(mat2str(gpm.delta_d_v_d_bw(1,:)))
disp(mat2str(gpm.delta_d_v_d_bw(2,:)))
disp(mat2str(gpm.delta_d_v_d_bw(3,:)))
disp('Numerical   ')
disp(mat2str(delta_d_v_d_bw_num(1,:)))
disp(mat2str(delta_d_v_d_bw_num(2,:)))
disp(mat2str(delta_d_v_d_bw_num(3,:)))
disp(' ')
disp('Jacobian d_p gyr bias')
disp('Analytical  ')
disp(mat2str(gpm.delta_d_p_d_bw(1,:)))
disp(mat2str(gpm.delta_d_p_d_bw(2,:)))
disp(mat2str(gpm.delta_d_p_d_bw(3,:)))
disp('Numerical   ')
disp(mat2str(delta_d_p_d_bw_num(1,:)))
disp(mat2str(delta_d_p_d_bw_num(2,:)))
disp(mat2str(delta_d_p_d_bw_num(3,:)))



% Acc bias Jacobian
gpm_dbf = cell(3,1);
delta_d_R_d_bf_num = zeros(3);
delta_d_v_d_bf_num = zeros(3);
delta_d_p_d_bf_num = zeros(3);
for i = 1:3
    imu_data_shift = imu_data;
    imu_data_shift.acc(:,i) = imu_data_shift.acc(:,i) + diff_delta;
    gpm_dbf{i} = Gpm(imu_data_shift.acc,...
            imu_data_shift.gyr,...
            imu_data_shift.time,...
            ground_truth.start_time,...
            ground_truth.end_time,...
            quantum,...
            simulation_opt.acc_std,...
            simulation_opt.gyr_std,...
            simulation_opt.one_axis);
    delta_d_v_d_bf_num(:,i) = (gpm_dbf{i}.d_v - gpm.d_v)/diff_delta;
    delta_d_p_d_bf_num(:,i) = (gpm_dbf{i}.d_p - gpm.d_p)/diff_delta;
end


disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('Jacobian d_v acc bias')
disp('Analytical  ')
disp(mat2str(gpm.delta_d_v_d_bf(1,:)))
disp(mat2str(gpm.delta_d_v_d_bf(2,:)))
disp(mat2str(gpm.delta_d_v_d_bf(3,:)))
disp('Numerical   ')
disp(mat2str(delta_d_v_d_bf_num(1,:)))
disp(mat2str(delta_d_v_d_bf_num(2,:)))
disp(mat2str(delta_d_v_d_bf_num(3,:)))
disp(' ')
disp('Jacobian d_p acc bias')
disp('Analytical  ')
disp(mat2str(gpm.delta_d_p_d_bf(1,:)))
disp(mat2str(gpm.delta_d_p_d_bf(2,:)))
disp(mat2str(gpm.delta_d_p_d_bf(3,:)))
disp('Numerical   ')
disp(mat2str(delta_d_p_d_bf_num(1,:)))
disp(mat2str(delta_d_p_d_bf_num(2,:)))
disp(mat2str(delta_d_p_d_bf_num(3,:)))
