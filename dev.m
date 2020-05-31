% Repository GPM - Gaussian Preintegrated Measurements
% This code is released under the MIT License.
% Copyright 2020 Cedric Le Gentil
%
% This code contains steps used for the implementation of the GPMs


%% SE kernel derivations

% Declare the symbolic variables
sf2 = sym('sf2',1,'real');
x = sym('x',1,'real');
z = sym('z',1,'real');
l2 = sym('l2',1,'real');
a = sym('a',1,'real');
b = sym('b',1,'real');
c = sym('c',1,'real');

% Kernel equation
k = sf2 * exp( - (x - z)*(x - z)/(2*l2) );


% Integration
k_int_left = int(k, x, [a b])
k_int_both = int(k_int_left, z, [a b]);
k_int_both = subs(k_int_both, z, x)

% Double integration
k_v = int(k,x,a,x);
k_v_int = int(k_v, x, [a b])
k_v_int_both = int(k_v_int,z, [a b]);
k_v_int_both = subs(k_v_int_both, z, x)

% For time-shift
k_int_left = int(k, x, [a a+c]);
k_int_left_dt = jacobian(k_int_left,a);
k_int_left_dt = simplify(k_int_left_dt)

k_v = int(k, x, a, x);
k_v_int = int(k_v, x, [a a+c]);
k_v_int_dt = jacobian(k_v_int,a);
k_v_int_dt = simplify(k_v_int_dt)

% For data time-shift
k_int_left = int(k, x, [a b]);
k_int_left_dt = jacobian(k_int_left,a);
k_int_left_dt = simplify(k_int_left_dt)
