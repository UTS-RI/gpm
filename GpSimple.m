% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
% 
% Function to infer a simple signal modelled with a Gaussian Process
function [varargout] = GpSimple(hyp, x, y, x_test)

% Read hyperparameters
sf2 = exp(2*hyp.cov(2));
l2 = exp(2*hyp.cov(1));

% Reshaping data for vector operations
xx = repmat(x, 1, length(x_test));
xx_test = repmat(x_test, 1,length(x))';
xxx = repmat(x, 1, length(x));

% Kernel first integral
ks = sf2 * exp( - (xx_test - xx).^2/(2*l2) );

% Value and associated variance nference
inverted_bit = inv(sf2 * exp( - (xxx - xxx').^2/(2*l2) ) + diag(ones(length(x),1)*exp(2*hyp.lik)));
my_alpha = inverted_bit * y;
mu = ks'*my_alpha;
cov_int = exp(2*hyp.lik)* sum((ks'*inverted_bit).^2,2);



varargout = {mu, cov_int};

