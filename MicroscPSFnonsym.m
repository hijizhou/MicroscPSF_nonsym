function PSF = MicroscPSFnonsym(params)
%MICROPSF Compute the 3D non-cylindrically symmetric PSF model described by Gibson and Lanni (JOSA 1992).
%   PSF = MICROPSF(params) return a 3D PSF given parameters
%
%   Parameters include:
%   (1) image properties
%           'size'  -  the size of the 3D PSF, e.g. params.size = [256 256 128];
%   (2) precision control
%           'numBasis'      - the number of approximation basis, default '100'
%           'numSamp'       - the number of sampling to determine the basis
%                             coefficients, default '1000'
%           'overSampling'  - the oversampling ratio, default 2
%   (3) microscope parameters
%        'NA'        - numerical aperture of the microscope, default 1.4
%        'lambda'    - Emission wavelength in vacuum, default 610nm
%        'M'         - magnification factor, default 100
%        'ns'        - specimen refractive index (RI), default 1.33
%        'ng0'       - coverslip RI, design value, default 1.5
%        'ng'        - coverslip RI, experimental, default 1.5
%        'ni0'       - immersion RI, design value, default 1.5
%        'ni'        - immersion RI, experimental, defualt 1.5
%        'ti0'       - working distance, design, default 150um
%        'tg0'       - coverslip thickness, design value, default 170um
%        'tg'        - coverslip thickness, experimental, default 170um
%        'resLateral' - lateral pixel size, default 100nm
%        'resAxial'  - axial pixel size, default 250nm
%        'pZ'        - position of particle, default 2000nm
%
%   Reference:
%       [1] Gibson, S.F. & Lanni, F., 1992.
%           Experimental test of an analytical model of aberration in an
%           oil-immersion objective lens used in three-dimensional light
%           microscopy. JOSA A, 9(1), pp.154-166.
%       [2] Li, J., Xue, F. and Blu, T. Fast and accurate 3D PSF
%           computation for fluorescence microscopy. JOSA A. Accepted.
%
%   See also AUX_BESSEL, AUX_SHOWPSF
%
%   Acknowledgement: PSFgenerator (http://bigwww.epfl.ch/algorithms/psfgenerator/)
%
%   Copyright © Jizhou Li, Feng Xue and Thierry Blu, 2018
%   Update date: 5 July, 2018
%   - add the non-cylindrically symmetric support
%   Note: the code is not fully tested. 

warning off;
if ~isfield(params,'size')
    error('Please set the size of PSF model');
end

size = params.size;
params.nx = size(1);
params.ny = size(2);
params.nz = size(3);

%% default parameters
if ~isfield(params,'numBasis')
    params.numBasis = 100;
end
if ~isfield(params,'numSamp')
    params.numSamp = 1000;
end
if ~isfield(params,'fastcom')
    params.fastcom = 0;
end
if ~isfield(params,'overSampling')
    params.overSampling = 2;
end
if ~isfield(params,'NA')
    params.NA = 1.4;
end
if ~isfield(params,'lambda')
    params.lambda = 610e-9;
end
if ~isfield(params,'M')
    params.M = 100;
end
if ~isfield(params,'ns')
    params.ns = 1.33;
end
if ~isfield(params,'ng0')
    params.ng0 = 1.5;
end
if ~isfield(params,'ng')
    params.ng = 1.5;
end
if ~isfield(params,'ni0')
    params.ni0 = 1.5;
end
if ~isfield(params,'ni')
    params.ni = 1.5;
end
if ~isfield(params,'ti0')
    params.ti0 = 150e-6;
end
if ~isfield(params,'tg0')
    params.tg0 = 170e-6;
end
if ~isfield(params,'tg')
    params.tg = 170e-6;
end
if ~isfield(params,'resLateral')
    params.resLateral = 200e-9;
end
if ~isfield(params,'resAxial')
    params.resAxial = 250e-9;
end
if ~isfield(params,'pZ')
    params.pZ = 2000e-9;
end

x0 = (params.nx-1)/2;
y0 = (params.ny-1)/2;
xp = x0; yp=y0;
maxRadius = round(sqrt((params.nx - x0).^2 + (params.ny - y0).^2)) + 1;
R = [0:params.overSampling*maxRadius-1]./params.overSampling;

a = 0;
b = min([1, params.ns/params.NA, params.ni/params.NA, params.ni0/params.NA,...
    params.ng0/params.NA,params.ng/params.NA]);

L = params.numSamp;
Rho = linspace(a,b,L)';

%% 1. approximate function exp(iW) as Bessel series
NN = params.numBasis;
k0 = 2*pi/params.lambda;
r = R*params.resLateral;
A = k0*params.NA*r;
A2 = A.^2;
Ab = A.*b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
an = (3*[1:NN]-2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

anRho = bsxfun(@times,an,Rho);
J = besselj(0, anRho);

J0A = besselj(0, Ab);
J1A = A.*besselj(1, Ab);

anJ0A = bsxfun(@times,J0A,an');
anb = an.*b;
an2 = an.^2;
B1anb = besselj(1, anb);
B0anb = besselj(0, anb);

Ele = bsxfun(@times,anJ0A,B1anb') - bsxfun(@times,J1A,B0anb');
domin = bsxfun(@minus, an2', A2);
Ele = Ele.*b./domin;

C1 = params.ns*params.pZ;
% C2 = params.ni*(Ti - params.ti0);
C3 = params.ng*(params.tg - params.tg0);
Z = params.resAxial*([0:params.nz-1] - ((params.nz - 1.0) / 2.0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Core
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalA = 181; % could be 361, finer
ANG = floor(linspace(-90, 90, totalA));
ANGindex = 1:totalA;
theta = ANG*pi/180;
[X,Y] = meshgrid(0:params.nx-1,0:params.ny-1);
rPixel = sqrt((X-xp).^2 + (Y-yp).^2);

theta2 = atand((Y-yp)./(X-xp));
theta2index = interp1(ANG,ANGindex, theta2);

for zi = 1: params.nz
    
    OPDs = C1*sqrt(1-(params.NA*Rho/params.ns).^2);
    % OPDi = bsxfun(@times, C2,sqrt(1-(params.NA*Rho/params.ni).^2));
    OPDg = C3*sqrt(1-(params.NA*Rho/params.ng).^2);
  
    deltaZ = params.pZ - Z(zi);
    
    OPD = bsxfun(@times, deltaZ*Rho.^2, cos(theta));
    
    % determine the coefficients
    W = k0*OPD;
    Ffun = cos(W) + 1i*sin(W);
    Ci = J\Ffun;
    
    %% 2. get PSF in each slice
    ciEle = Ele'*Ci;
    PSF0 = abs(ciEle).^2;
    
    %% 3. obtain the 2D slice
    slice = interp2(PSF0, theta2index, rPixel+1);
%     imshow(slice,[]); drawnow;
    PSF(:,:,zi) = slice;
end

PSF = PSF./max(PSF(:));

end
