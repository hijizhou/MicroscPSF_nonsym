%   Update date: 5 July, 2018

addpath('Utilities/');

clear; clc;
params.size = [128 128 64];
% params.NA = 1;
params.pZ = 0;
params.ns = 1.5;

tic;
PSF = MicroscPSFnonsym(params);
t = toc;

% icy_im3show([PSF]);
disp(['Running time = ' num2str(t) 's']);

