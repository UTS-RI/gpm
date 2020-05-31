% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% Function to generate a simulated trajectory and the associated IMU
% measurements
function [ ground_truth, imu_data ] = SimulateTraj( options )


    % Generate the gravity direction (elevation and azymuth)
    if options.one_axis
        g_simu = [0, 0];
        disp(['Simulate a 1-axis rotation trajectory with the '...
            options.traj_profile ' profile']);
    else
        g_simu = randn(1,2);
        disp(['Simulate a 3-axis rotation trajectory with the '...
            options.traj_profile ' profile']);
    end

    % Create gravity vector from elevation and azymuth
    g = options.gravity_magnitude*[sin(g_simu(1))*cos(g_simu(2));...
                                   sin(g_simu(1))*sin(g_simu(2));...
                                   cos(g_simu(1))];
    
                               
                               
    % Define the trajectory parameters depending on the motion type
    if strcmp(options.traj_profile, 'slow')
        vel_amplitude_min = 0.3;
        vel_amplitude_max = 2;
        vel_freq_min = 0.05;
        vel_freq_max = 0.1;
        rot_amplitude_min = 0.5;
        rot_amplitude_max = 0.75;
        rot_freq_min = 0.05;
        rot_freq_max = 0.3;
    elseif strcmp(options.traj_profile, 'fast')
        vel_amplitude_min = 3;
        vel_amplitude_max = 5;
        vel_freq_min = 0.7;
        vel_freq_max = 1.5;
        rot_amplitude_min = 0.75;
        rot_amplitude_max = 1.5;
        rot_freq_min = 0.25;
        rot_freq_max = 0.75;
    else
        disp('ERROR: The trajectory profile is not known');
    end

    
    % Create the sine functions characteristics from the previous
    % parameters
    nb_sine = 3;
    vel_sine = cell(nb_sine,3);
    for i = 1:nb_sine
        v_f_increment = (vel_freq_max - vel_freq_min)/nb_sine;
        v_f_min = vel_freq_min + (i-1)*v_f_increment;
        v_f_max = vel_freq_min + (i)*v_f_increment;
        
        for j = 1:3
            vel_sine{i,j}.freq = (v_f_max - v_f_min)*rand() + v_f_min;
            vel_sine{i,j}.amp =...
                (vel_amplitude_max - vel_amplitude_min)*rand()...
                + vel_amplitude_min;
            
        end

    end
    
    rot_sine = cell(1,3);
    for j = 1:3
        rot_sine{j}.freq = (rot_freq_max - rot_freq_min)*rand()...
            + rot_freq_min;
        rot_sine{j}.amp =...
            (rot_amplitude_max - rot_amplitude_min)*rand()...
            + rot_amplitude_min;        
    end
    if options.one_axis
        rot_sine{2}.amp = 0;
        rot_sine{3}.amp = 0;
    end
        
    % Generate the timeline
    start_time = 100*rand();
    end_time = start_time + options.duration;
    imu_period = 1.0/options.imu_frequency;
    overlap_low = start_time - options.data_overlap;
    overlap_high = end_time + options.data_overlap;
    imu_start_time = floor(overlap_low/imu_period)*imu_period;
    imu_end_time = ceil(overlap_high/imu_period)*imu_period;
    imu_time = imu_start_time:imu_period:imu_end_time;
    dataset_size = length(imu_time);
    imu_time = [imu_time start_time end_time]';

    % Generate linear pos, vel and acc
    vel = zeros(dataset_size+2,3);
    acc = zeros(dataset_size+2,3);
    pos = zeros(dataset_size+2,3);
    % First deal with the x axis velocity offset
    vel(:,1) = options.x_vel_offset;
    pos(:,1) = options.x_vel_offset*imu_time;
    % Loop through the different sine functions
    for i = 1:nb_sine
        for j = 1:3
            vel(:,j) = vel(:,j) + vel_sine{i,j}.amp*cos(...
                2*pi*vel_sine{i,j}.freq*imu_time );
            acc(:,j) = acc(:,j) - vel_sine{i,j}.amp*sin(...
                2*pi*vel_sine{i,j}.freq*imu_time )...
                * ( 2*pi*vel_sine{i,j}.freq );
            pos(:,j) = pos(:,j) + vel_sine{i,j}.amp*sin(...
                2*pi*vel_sine{i,j}.freq*imu_time)...
                /(2*pi*vel_sine{i,j}.freq);
        end
    end
    
   
    % Generate orientation and angular vel
    euler = ...
            [ rot_sine{1}.amp*sin(2*pi*rot_sine{1}.freq*imu_time),...
              rot_sine{2}.amp*sin(2*pi*rot_sine{2}.freq*imu_time),...
              rot_sine{3}.amp*sin(2*pi*rot_sine{3}.freq*imu_time)...
            ];

    ang_vel = ...
        [                                                                                   2*rot_sine{3}.amp*rot_sine{3}.freq*pi*cos(2*pi*rot_sine{3}.freq*imu_time) - 2*rot_sine{1}.amp*rot_sine{1}.freq*pi*sin(rot_sine{2}.amp*sin(2*pi*rot_sine{2}.freq*imu_time)).*cos(2*pi*rot_sine{1}.freq*imu_time),...
         2*rot_sine{2}.amp*rot_sine{2}.freq*pi*cos(2*pi*rot_sine{2}.freq*imu_time).*cos(rot_sine{3}.amp*sin(2*pi*rot_sine{3}.freq*imu_time)) + 2*rot_sine{1}.amp*rot_sine{1}.freq*pi*sin(rot_sine{3}.amp*sin(2*pi*rot_sine{3}.freq*imu_time)).*cos(2*pi*rot_sine{1}.freq*imu_time).*cos(rot_sine{2}.amp*sin(2*pi*rot_sine{2}.freq*imu_time)),...
         2*rot_sine{1}.amp*rot_sine{1}.freq*pi*cos(2*pi*rot_sine{1}.freq*imu_time).*cos(rot_sine{2}.amp*sin(2*pi*rot_sine{2}.freq*imu_time)).*cos(rot_sine{3}.amp*sin(2*pi*rot_sine{3}.freq*imu_time)) - 2*rot_sine{2}.amp*rot_sine{2}.freq*pi*sin(rot_sine{3}.amp*sin(2*pi*rot_sine{3}.freq*imu_time)).*cos(2*pi*rot_sine{2}.freq*imu_time)...
        ];

    % Create accelerometer measurements from acceleration, gravity and
    % orientation
    imu_acc = zeros(dataset_size,3);
    imu_rot = zeros(3,3,dataset_size);
    for i = 1 : dataset_size
        rot = eul2rotm(euler(i,:));
        imu_acc(i,:) = (rot'*(acc(i,:)'))' - ( (rot')*g)';
        imu_rot(:,:,i) = rot;
    end
    imu_data.time = imu_time(1:(end-2));
    imu_data.acc = imu_acc + options.acc_std*randn(size(imu_acc));
    imu_data.gyr = ang_vel(1:(end-2),:)...
        + options.gyr_std*randn(size(imu_acc));

    if options.one_axis
        imu_data.gyr(:, 1:2) = 0;
    end
    
        
    % Store the ground truth trajectories
    ground_truth.vel_start = vel(end-1,:)';
    ground_truth.vel_end = vel(end,:)';
    ground_truth.pos_start = pos(end-1,:)';
    ground_truth.pos_end = pos(end,:)';
    ground_truth.rot_start = eul2rotm(euler(end-1,:));
    ground_truth.rot_end = eul2rotm(euler(end,:));
    d_t = end_time - start_time;
    ground_truth.imu_time = imu_time(1:(end-2));
    ground_truth.imu_pos = pos(1:(end-2),:);
    ground_truth.imu_vel = vel(1:(end-2),:);
    ground_truth.imu_rot = imu_rot;
    ground_truth.start_time = start_time;
    ground_truth.end_time = end_time;
    ground_truth.gravity = g;

    % Visualisation of the trajectory if wanted
    if options.visualisation
        figure(1010)
        plot3(pos(1:(end-2),1), pos(1:(end-2),2), pos(1:(end-2),3))
        hold on
        plot3(ground_truth.pos_start(1), ground_truth.pos_start(2), ground_truth.pos_start(3), 'ob')
        plot3(ground_truth.pos_end(1), ground_truth.pos_end(2), ground_truth.pos_end(3), 'or')
        hold off
        legend('Trajectory', 'Start pos', 'End pos');
        title('Simulated trajectory visualisation');
        axis('equal')
    end
    
    % Compute and store the ground truth preintegrated measurements
    ground_truth.d_R = ground_truth.rot_start' * ground_truth.rot_end;
    ground_truth.d_v = ground_truth.rot_start'...
        * (ground_truth.vel_end - ground_truth.vel_start - (d_t*g) );
    ground_truth.d_p = ground_truth.rot_start'...
        * (ground_truth.pos_end - ground_truth.pos_start ...
            - (d_t*ground_truth.vel_start) - (d_t*d_t*g/2));


end
