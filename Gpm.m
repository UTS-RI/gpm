% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% Computation of the GPMs and the associated postintegration Jacobians
function gpm_out = Gpm( acc, gyr, time, start_time,  query_time,...
    quantum, acc_sd, gyr_sd, one_axis_flag)


    % Parameters for the GP fitting
    likelihood = 'likGauss';
    inference = 'infGaussLik';
    cov_function = {@covSEiso};
    mean_function = @meanZero;
    ncg = 50;
    
    % Prepare output container
    gpm_out = cell(length(query_time),1);
    for i = 1:length(query_time)
        gpm_out{i}.d_R = eye(3);
        gpm_out{i}.d_v = zeros(3,1);
        gpm_out{i}.d_p = zeros(3,1);
        gpm_out{i}.delta_d_R_d_t = zeros(3,1);
        gpm_out{i}.delta_d_v_d_t = zeros(3,1);
        gpm_out{i}.delta_d_p_d_t = zeros(3,1);
        gpm_out{i}.delta_d_v_d_bf = zeros(3,3);
        gpm_out{i}.delta_d_p_d_bf = zeros(3,3);
        gpm_out{i}.delta_d_R_d_bw = zeros(3,3);
        gpm_out{i}.delta_d_v_d_bw = zeros(3,3);
        gpm_out{i}.delta_d_p_d_bw = zeros(3,3);
        gpm_out{i}.cov = zeros(9,9);
    end
    data_d_t = zeros(size(acc));
    rot_d_t = zeros(size(acc));
    data_d_bw = cell(3,1);
    rot_d_bw = cell(3,1);
    data_d_bf = cell(3,1);
    for i = 1:3
        data_d_bw{i} = zeros(size(acc));
        rot_d_bw{i} = zeros(size(acc));
        data_d_bf{i} = zeros(size(acc));
    end


    % Organise the gyr data to be able to loop between the sensors axis
    ytr = [...
        gyr(:,1),...
        gyr(:,2),...
        gyr(:,3)];
        
    xtr = time;
    numerical_delta = 0.0001;
    
    
    % Prepare data structure for the 1-axis rotation case
    one_d_time = [time; start_time; query_time];
    d_r = zeros(length(one_d_time), 3);
    d_r_cov = zeros(length(one_d_time), 3);
    
    
    
    
    % To store the signals' mean and learnt hyper-parameters for the 3-axis
    % rotations inference further
    m = zeros(3,1);
    hyp_save = cell(3,1);
    
    
    % Fit GP to gyr data and infer integral at acc timestamps and query
    % time if 1-axis rotations
    for i = 1:3
               

        % Centre the data around zero
        m(i) = sum(ytr(:,i))/length(ytr(:,i));
        ytr(:,i) = ytr(:,i) - m(i);

        % Train hyper-parameter
        hyp = TrainHyp(ncg, inference, mean_function, cov_function,...
                            likelihood, xtr, ytr(:,i), gyr_sd);
        hyp_save{i} = hyp;

        
                
        % If 1-axis rotations, directly infer the preintegration rotational
        % measurements
        if one_axis_flag

            if i == 3
                [d_r(:,i), d_r_cov(:,i),...
                delta_d_r_dbw, delta_d_r_dt,...
                delta_data_d_r_dt] = GpIntegral2(hyp,...
                    xtr, ytr(:,i), start_time, one_d_time,...
                    ones(size(ytr(:,i))) );
                d_r(:,i) = d_r(:,i) + (m(i)*(one_d_time - start_time));
                delta_data_d_r_dt = delta_data_d_r_dt - m(i);
            else
                d_r(:,i) = zeros(size(d_r(:,i)));
                d_r_cov(:,i) = zeros(size(d_r_cov(:,i)));
                delta_d_r_dbw = zeros(size(one_d_time));
                delta_d_r_dt = zeros(size(one_d_time));
                delta_data_d_r_dt = zeros(size(one_d_time));
            end

            % Extract the Jacobians for the output
            for j = 1:length(query_time)
                index = length(one_d_time) - length(query_time) + j;
                gpm_out{j}.delta_d_R_d_t(i) = delta_d_r_dt(index);
                gpm_out{j}.delta_d_R_d_bw(i,i) = delta_d_r_dbw(index);
            end
            % Get the Jacobian bits of the data preintegration measurements
            for j = 1:length(time)
                rot_d_t(j,i) = delta_data_d_r_dt(j);
                rot_d_bw{i}(j,i) = delta_d_r_dbw(j);
            end
        end
    end
    %%%

    
    % If three axis rotation, orientation computed numerically
    % and rotate the acceleration measurements
    acc_rot = zeros(length(time), 3);
    acc_rot_d_t = zeros(length(time), 3);
    acc_rot_d_bw = cell(3,1);
    for i = 1:3
        acc_rot_d_bw{i} = zeros(length(time), 3);
    end
    if ~one_axis_flag
        % Interpolate the gyroscope measurements for the IMU data
        % time-stamp, the query points, and the time-stamps needed for
        % time-shift numerical Jacobians
        t = min(time):quantum:max(time);
        t = [t start_time query_time' time'];
        [t, id_sorting] = sort(t);
        t = t';
        t_extended = [t;...
            (start_time+numerical_delta);...
            (query_time+numerical_delta)];
        
        gyr_dense = zeros(length(t_extended),3);
        gyr_dense_cov = zeros(length(t_extended),3);
        for i = 1:3
            [gyr_dense(:,i), gyr_dense_cov(:,i)] = GpSimple(...
                hyp_save{i}, xtr, ytr(:,i), t_extended);
            gyr_dense(:,i) = gyr_dense(:,i) + m(i);
        end
        gyr_dense_d_t(id_sorting,:) = gyr_dense(1:length(t),:);
        t_dt(id_sorting) = t;
        t_dt(end-length(time)-length(query_time))...
            = (start_time+numerical_delta);   
        t_dt((end-length(time)-length(query_time)+1):...
            end-length(time))...
            = (query_time+numerical_delta);
        [t_dt, id_sorting_dt] = sort(t_dt);
        gyr_dense_d_t = gyr_dense_d_t(id_sorting_dt, :);
        gyr_dense = gyr_dense(1:length(t),:);
        gyr_dense_cov = gyr_dense_cov(1:length(t),:);
        
        
        % Preintegrate
        [d_R_raw, d_R_cov] = RotationNumericalPreintegration(...
                t, gyr_dense, gyr_dense_cov, start_time,query_time(end));
        
            
        % Preintegrations for the numerical Jacobians
        d_R_raw_dt = RotationNumericalPreintegration(...
            t_dt, gyr_dense_d_t);
        
        d_R_raw_d_bw = cell(3,1);
        for i = 1:3
            gyr_temp = gyr_dense;
            gyr_temp(:,i) = gyr_temp(:,i) + numerical_delta;
            d_R_raw_d_bw{i} = RotationNumericalPreintegration(...
                t, gyr_temp);
        end
            
        % Demux the rotations to later project the acc measurements    
        d_R = zeros(3,3,length(t));
        d_R(:,:,id_sorting) = d_R_raw;
        d_R_acc = d_R(:,:,end-(length(time)-1):end);
        d_R_start = d_R(:,:,end-length(time)-length(query_time));
                    
        d_R_d_t = zeros(3,3,length(t));
        d_R_d_t(:,:,id_sorting_dt) = d_R_raw_dt;
        d_R_acc_d_t = d_R_d_t(:,:,end-(length(time)-1):end);
        d_R_start_d_t = d_R_d_t(:,:,end-length(time)-length(query_time));
        
        d_R_d_bw = cell(3,1);
        d_R_acc_d_bw = cell(3,1);
        d_R_start_d_bw = cell(3,1);
        for i = 1:3
            d_R_d_bw{i} = zeros(3,3,length(t));
            d_R_d_bw{i}(:,:,id_sorting) = d_R_raw_d_bw{i};
            d_R_acc_d_bw{i} = d_R_d_bw{i}(:,:,end-(length(time)-1):end);
            d_R_start_d_bw{i} = d_R_d_bw{i}(...
                :,:,end-length(time)-length(query_time));
        end
        
        
        % Fill the output data structure with the rotational part
        for i = 1:length(query_time)
            gpm_out{i}.d_R = d_R_start'...
                *d_R(:,:,end-length(time)-length(query_time) + i);
            gpm_out{i}.cov(1:3,1:3) =...
                d_R_cov(:,:,end-length(time)-length(query_time) + i);
            gpm_out{i}.delta_d_R_d_t =...
                LogMap(gpm_out{i}.d_R'*d_R_start_d_t'...
                *d_R_d_t(:,:,end-length(time)-length(query_time) + i))...
                /numerical_delta;
            for j = 1:3
                data_d_bw_temp = LogMap(...
                    gpm_out{i}.d_R'*d_R_start_d_bw{j}'...
                    *d_R_d_bw{j}(:,:,end-length(time)-length(query_time)...
                        + i))...
                    /numerical_delta;
                gpm_out{i}.delta_d_R_d_bw(:,j) = data_d_bw_temp;
            end
        end
        
        % Project the acc data and get the associated Jacobians
        for i = 1:length(time)
            temp_rot = d_R_start' * d_R_acc(:,:,i);
            acc_rot(i,:) = (temp_rot * acc(i,:)')';
            acc_rot_d_t(i,:) = (d_R_start_d_t'...
                * d_R_acc_d_t(:,:,i) * acc(i,:)')';
            for j = 1:3
                acc_rot_d_bw{j}(i,:) = (d_R_start_d_bw{j}'...
                    * d_R_acc_d_bw{j}(:,:,i) * acc(i,:)')';
                data_d_bw_temp = (acc_rot_d_bw{j}(i,:) - acc_rot(i,:))...
                /numerical_delta;
                data_d_bw{1}(i,j) = data_d_bw_temp(1);
                data_d_bw{2}(i,j) = data_d_bw_temp(2);
                data_d_bw{3}(i,j) = data_d_bw_temp(3);
                
                data_d_bf{j}(i,:) = temp_rot(j,:);
            end
            data_d_t(i,:) = (acc_rot_d_t(i,:) - acc_rot(i,:))...
                /numerical_delta;
        end
        
        
    % Else is the rotation is only around 1 axis exponentiate the previous
    % GP output to obtain the rotation matrices and rotation accelerometers
    % measurements
    else
        % Project the acc data and associated Jacobians
        R_start = ExpMap(d_r(end-length(query_time),:));
        for i = 1:length(time)
            rot_temp = (R_start' * ExpMap(d_r(i,:)));
            acc_rot(i,:) = (rot_temp * acc(i,:)')';
            for j = 1:3
                data_d_bf{j}(i,:) = rot_temp(j,:);
                q = zeros(3,1);
                q(j) = 1;
                temp_1 = rot_temp'*q*acc(i,:);
                temp_2 = ExpMapJacobian(zeros(3,1));
                temp_3 = [rot_d_bw{1}(i,:);...
                            rot_d_bw{2}(i,:);...
                            rot_d_bw{3}(i,:)];
                temp_4 = rot_d_t(i,:)';
                data_d_bw{j}(i,:) = temp_1(:)' * temp_2 * temp_3;
                data_d_t(i,j) = temp_1(:)' * temp_2 * temp_4;
            end
        end
        
        
        % Fill the output data structure with the rotational part
        for i = 1:length(query_time)
            gpm_out{i}.d_R = R_start'...
                *ExpMap(d_r(end-length(query_time)+i,:));
            gpm_out{i}.cov(1:3,1:3) =...
                diag(d_r_cov(end-length(query_time)+i,:));
        end
        
    end
    
    
    
 

    % Orgarnise the acc data to loop through
    ytr = [...
        acc_rot(:,1),...
        acc_rot(:,2),...
        acc_rot(:,3)];


    
    
    
    % Create the velocity and position preintegrated measurements
    for i = 1:3


        % Centre the data around zero
        m(i) = sum(ytr(:,i))/length(ytr(:,i));
        ytr(:,i) = ytr(:,i) - m(i);

        % Train hyper-parameter
        hyp = TrainHyp(ncg, inference, mean_function, cov_function,...
                            likelihood, xtr, ytr(:,i), acc_sd);
        
        
        
        
        
        % Integrals inferences
        [d_v, d_v_cov, d_p, d_p_cov,...
            delta_d_v_d_bw, delta_d_p_d_bw,...
            delta_d_v_d_bf, delta_d_p_d_bf,...
            delta_d_v_d_t, delta_d_p_d_t] = GpIntegral2(hyp,...
                        xtr, ytr(:,i), start_time, query_time,...
                        data_d_bw{i}, data_d_bf{i},...
                        data_d_t(:,i));
        
        
        % Fill up the output data structure
        for j = 1:length(query_time)
            gpm_out{j}.d_v(i) = d_v(j)...
                + (m(i)*(query_time(j) - start_time));
            gpm_out{j}.d_p(i) = d_p(j)...
                + (m(i)*((query_time(j) - start_time).^2)/2);
            gpm_out{j}.cov(3+i,3+i) = d_v_cov(j);
            gpm_out{j}.cov(6+i,6+i) = d_p_cov(j);
            gpm_out{j}.delta_d_p_d_t(i) = delta_d_p_d_t(j);
            gpm_out{j}.delta_d_v_d_t(i) = delta_d_v_d_t(j);
            gpm_out{j}.delta_d_p_d_bf(i,:) = delta_d_p_d_bf(j,:);
            gpm_out{j}.delta_d_v_d_bf(i,:) = delta_d_v_d_bf(j,:);
            gpm_out{j}.delta_d_p_d_bw(i,:) = delta_d_p_d_bw(j,:);
            gpm_out{j}.delta_d_v_d_bw(i,:) = delta_d_v_d_bw(j,:);
        end
        
        
        
    end    
    
    
    
    % Reshape the output if only one element
    if length(query_time) == 1
        gpm_out = gpm_out{1};
    end
    
end


function hyp = TrainHyp(ncg, inference, mean_function, cov_function,...
                            likelihood, xtr, ytr, obs_std)


        % Prior on the hyper-parameters
        sf = std(ytr);
        ell = 5*(xtr(end) - xtr(1) ) / (length(xtr) - 1);
        % If the signal standard deviation is smaller than sensor noise
        % do not train hyper-parameters (signal approximatively constant,
        % use prior)
        max_nb_loops = 50;
        if sf < 1.1*obs_std
            sf = 0.2*obs_std;
            max_nb_loops = 0;
        end
        hyp0.cov  = log([ell;sf]);
        hyp0.lik  = log(obs_std);
        hyp = hyp0;
        
        

        
        % Learn hyper-parameters (there is a simple mechanism to check if
        % the hyper-parameters make sense)
        % "max_nb_loops" corresponds the maximum number of tries to learn
        % the hyper-parameters. The method implemented is not bullet-proof
        % but seems to work well enough :))
        if(max_nb_loops ~= 0)
            hyp = minimize(hyp0,'gp', -ncg, inference, mean_function,...
                cov_function, likelihood, xtr, ytr);
            
            if(max_nb_loops > 1)
                counter = 1;
                loop = max_nb_loops > 1;
                ins = NaN(max_nb_loops,1);
                hyps = cell(max_nb_loops,1);
                while loop
                    [test_out, test_cov] = gp(hyp, inference,...
                        mean_function, cov_function,...
                        likelihood, xtr, ytr, xtr);

                    low_1 = test_out - sqrt(test_cov);
                    hi_1 = test_out + sqrt(test_cov);
                    low_5 = test_out - 5*sqrt(test_cov);
                    hi_5 = test_out + 5*sqrt(test_cov);

                    in_1 = sum( (ytr > low_1) & (ytr < hi_1) )/length(xtr);
                    in_5 = sum( (ytr > low_5) & (ytr < hi_5) )/length(xtr);
                    ins(counter) = in_1;
                    hyps{counter} = hyp;
                    in_5_thr = (length(xtr) - 2)/length(xtr);

                    if (in_1 < 0.95) && (in_1 > 0.5) && (in_5 > in_5_thr) 
                        loop = false;
                    elseif counter >= max_nb_loops
                        loop = false;
                        ins_error = abs(ins - 0.6);

                        [min_ins, id_min] = min(ins_error);
                        hyp = hyps{id_min};
                        disp(['Ins = ' num2str(ins(id_min))]);
                    end

                    if loop == true
                        new_ell = max(exp(hyp.cov(1))...
                            + 0.1*randn(), 0.00001);
                        new_sf = max(exp(hyp.cov(2))...
                            + 0.1*randn(), 0.00001);
                        hyp0.cov = log([new_ell; new_sf]);
                        hyp = minimize(hyp0,'gp', -ncg,...
                            inference, mean_function, cov_function,...
                            likelihood, xtr, ytr);
                    end
                    counter = counter + 1;
                
                end
            end
        end
end

% Preintegration function for rotational part
function [varargout] = RotationNumericalPreintegration(...
                                    t, gyr, gyr_var, start_time, end_time)

    d_t = t(2:end) - t(1:(end-1));
    if nargout > 1
        d_R_cov = zeros(3,3,length(t));
    end
    d_R_raw = zeros(3,3,length(t));
    d_R_raw(:,:,1) = eye(3);
    d_R = eye(3);
    cov_pre_int = zeros(3);
    for i = 1:(length(t)-1)
        if d_t(i) ~= 0
            e_R = ExpMap((gyr(i,:)) * d_t(i));
            if nargout > 1
                if (t(i)>= start_time) && (t(i)< end_time)
                    cov_imu = diag(gyr_var(i,:));

                    j_r = So3RightJacobian((gyr(i,:)) * d_t(i));

                    A = e_R';
                    B = j_r*d_t(i);
                    cov_pre_int...
                        = (A * cov_pre_int * A') + (B * cov_imu * B');
                end
                d_R_cov(:,:,i+1) = cov_pre_int;
            end

            d_R = d_R * e_R;
            d_R_raw(:,:,i+1) = d_R;
        else
            d_R_raw(:,:,i+1) = d_R_raw(:,:,i);
            if nargout > 1
                d_R_cov(:,:,i+1) = d_R_cov(:,:,i);
            end
        end
    end

    if nargout > 1
        varargout = {d_R_raw, d_R_cov};
    else
        varargout = {d_R_raw};
    end
    
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
