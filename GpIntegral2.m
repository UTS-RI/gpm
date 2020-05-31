% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% Function to infer the integral and second integral of a signal modelled
% with a Gaussian Process
function [varargout] = GpIntegral2(hyp, x, y, a ,b,...
                                delta_y_d_bw, delta_y_d_bf , delta_y_d_t)

% Read hyperparameters
sf2 = exp(2*hyp.cov(2));
l2 = exp(2*hyp.cov(1));
sq_pi = pi^(1/2);

% Reshaping data for vector operations
xx = repmat(x, 1, length(b));
bb = repmat(b', length(x), 1);
xxx = repmat(x, 1, length(x));

% Kernel first integral
ks = -(2^(1/2)*sf2*sq_pi*(erf((2^(1/2)*(a - xx)*(1/l2)^(1/2))/2) - erf((2^(1/2)*(bb - xx)*(1/l2)^(1/2))/2)))/(2*(1/(l2))^(1/2));

% First integral and associated variance inference
inverted_bit = inv(sf2 * exp( - (xxx - xxx').^2/(2*l2) ) + diag(ones(length(x),1)*exp(2*hyp.lik)));
my_alpha = inverted_bit * y;
mu_int = ks'*my_alpha;
cov_int = exp(2*hyp.lik)* sum((ks'*inverted_bit).^2,2);


if ~(((nargout == 3) || (nargout == 4) || (nargout == 5)) && (nargin > 5))
    % Kernel second integral
    ks2 = (2^(1/2)*sf2*sq_pi*erf((2^(1/2)*(a - xx)*(1/l2)^(1/2))/2).*(a - bb)*(l2)^(1/2))/2 - (2^(1/2)*sf2*sq_pi*(l2)^(1/2)*(erf((2^(1/2)*(a - xx)*(1/l2)^(1/2))/2).*(a - xx) - erf((2^(1/2)*(bb - xx)*(1/l2)^(1/2))/2).*(bb - xx) + (2^(1/2)*exp(-(a - xx).^2/(2*l2))*(l2)^(1/2))/sq_pi - (2^(1/2)*exp(-(bb - xx).^2/(2*l2))*(l2)^(1/2))/sq_pi))/2;

    % Second integral and associated variance nference
    mu_int_2 = ks2'*my_alpha;
    cov_int_2 = exp(2*hyp.lik)* sum((ks2'*inverted_bit).^2,2);

    varargout{3} = mu_int_2;
    if nargout > 3
        varargout{4} = cov_int_2;
    end
end 





% Gyr bias jacobians
if nargin > 5
    mu_int_dbw = ks'*inverted_bit*delta_y_d_bw;
    if (nargout == 3) || (nargout == 4) || (nargout == 5)
        varargout{3} = mu_int_dbw;
    else
        mu_int_2_dbw = ks2'*inverted_bit*delta_y_d_bw;
        varargout{5} = mu_int_dbw;
        varargout{6} = mu_int_2_dbw;
    end
end


% Acc bias jacobians
if nargin > 6
    mu_int_dbf = ks'*inverted_bit*delta_y_d_bf;
    mu_int_2_dbf = ks2'*inverted_bit*delta_y_d_bf;
    varargout{7} = mu_int_dbf;
    varargout{8} = mu_int_2_dbf;
end


% Time-shift jacobians
if (nargin > 7) || ( (nargout > 3) && (nargin > 5) )
    ks_dt = sf2*exp(-(bb - xx).^2/(2*l2)) - sf2*exp(-(a - xx).^2/(2*l2));

    if nargin <= 7
        mu_int_dt = ks_dt'*my_alpha;
        varargout{4} = mu_int_dt;
    else
        mu_int_dt = ks_dt'*my_alpha + ks'*inverted_bit*delta_y_d_t;
        ks2_dt = -(2^(1/2)*sf2*sq_pi*erf((2^(1/2)*(a - xx)*(1/l2)^(1/2))/2) - 2^(1/2)*sf2*sq_pi*erf((2^(1/2)*(1/l2)^(1/2)*(bb - xx))/2) + 2*(b - a)*sf2*exp(-(a - xx).^2/(2*l2))*(1/l2)^(1/2))/(2*(1/l2)^(1/2));

        mu_int_2_dt = ks2_dt'*my_alpha + ks2'*inverted_bit*delta_y_d_t;
        varargout{9} = mu_int_dt;
        varargout{10} = mu_int_2_dt;
    end
end

% Data time-shift jacobian
if nargout == 5
    ks_dt = -sf2*exp(-(a - xx).^2/(2*l2));
    mu_int_dt = ks_dt'*my_alpha;
    varargout{5} = mu_int_dt;
end


varargout{1} = mu_int;
varargout{2} = cov_int;


