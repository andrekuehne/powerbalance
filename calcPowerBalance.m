function powerBalance = calcPowerBalance(e, h, sigma, GridVectors, uMatrix, iMatrix, sMatrix, forwardPower, varargin)
%calcPowerBalance: Calculate the power balance matrices of a multiport EM FDTD simulation with nPorts ports.
%   If used for a publication, please cite Kuehne et al. Power Balance and Loss Mechanism Analysis for RF Transmit Coil Arrays
%   https://pubmed.ncbi.nlm.nih.gov/25324179/
%   calcPowerBalance(e, h, sigma, GridVectors, uMatrix, iMatrix, sMatrix, forwardPower)
%   returns a struct containing all power balance matrices.
%
%   Input arguments:
%       e               - An array containing the 3D electric (E) fields of
%                         size [nPorts nx ny nz 3].
%       h               - An array containing the 3D magnetic (H) fields of
%                         size [nPorts nx ny nz 3].
%       sigma           - An array containing electric material conductivity of
%                         size [nx ny nz 3].
%       GridVectors     - A cell array containing the Yee cell grid vectors
%                         (Yee cell corner locations along each axis).
%                         Ex.: {xLine, yLine, zLine}, with xLine, yLine and
%                         zLine being of size [nx+1], [ny+1] and [nz+1],
%                         respectively (also see diagram below).
%                       
%       uMatrix         - An nLumpedElements x nPorts matrix containing the 
%                         lumped element voltages.
%       iMatrix         - An nLumpedElements x nPorts matrix containing the 
%                         lumped element currents.
%       sMatrix         - The scattering matrix of the system.
%       forwardPower    - The forward power used for the EM simulation.
%
%   Optional input arguments:
%       'farFieldOffset'  - Distance from the outer field boundaries for
%                           Poynting flux integration (default is 4).
%       'sigmaMasks'      - Cell array of binary masks for material loss.
%                           matrix discrimination. Each binary mask needs
%                           to be of size [nx ny nz 3].
%       'lumpedMasks'     - Cell array of binary masks for lumped element
%                           loss discrimination. Each binary mask needs
%                           to be of length nLumpedElements.
%
%   Arrangement of data arrays
%
%   The 3D field Data is expected to be distributed on a 3D Yee grid.
%   Yee cell coordinates need to be located at the cell vertices. For
%   the relative positioning of cell coordintes and field components
%   refer to the 2D Grid example:
%
%    (x1,y3)           (x2,y3)           (x3,y3)
%       |                 |
%       |                 |
%       |                 |
%   Ey(1,2)           Ey(2,2)
%       |                 |
%       |                 |
%       |                 |
%       |_ _ _Ex(1,2)_ _ _|_ _ _Ex(2,2)_ _ _
%    (x1,y2)           (x2,y2)           (x3,y2)
%       |                 |
%       |                 |
%       |                 |
%   Ey(1,1)           Ey(2,1)
%       |                 |
%       |                 |
%       |                 |
%       |_ _ _Ex(1,1)_ _ _|_ _ _Ex(2,1)_ _ _
%    (x1,y1)           (x2,y1)           (x3,y1)
%
%   The GridVectors property would in this case be
%   GridVectors = {[x1, x2, x3], [y1, y2, y3]};

%% Input parsing and sanity checks
narginchk(8,14);
varInput = varargin;
[nPorts, farFieldOffset, sigmaMasks, lumpedMasks] = parseInput(e, h, sigma, GridVectors, uMatrix, iMatrix, sMatrix, forwardPower, varInput);

%% Calculate matrices
%% forward
powerBalance.fwd = calcForwardPowerMatrix(nPorts, forwardPower);

%% coupling
powerBalance.cpl = calcCouplingLossMatrix(sMatrix, forwardPower);

%% lumped elements
powerBalance.lmp = calcLumpedLossMatrix(uMatrix, iMatrix);

if ~isempty(lumpedMasks)
   powerBalance.lmpMasked=cell(length(lumpedMasks),1);
   
   for k=1:length(lumpedMasks)
      currMask = (lumpedMasks{k}(:));
      currU = bsxfun(@times, uMatrix, currMask);
      currI = bsxfun(@times, iMatrix, currMask);
      powerBalance.lmpMasked{k} = calcLumpedLossMatrix(currU, currI);
       
   end
    
end

%% material loss
powerBalance.mat = calcMaterialLossMatrix(e, sigma, GridVectors{1}, GridVectors{2}, GridVectors{3});

% if there are material masks loop through them
if ~isempty(sigmaMasks)
    powerBalance.matMasked=cell(length(sigmaMasks),1);
    
    for k=1:length(sigmaMasks)
        
        powerBalance.matMasked{k} = calcMaterialLossMatrix(e, sigma.*sigmaMasks{k}, GridVectors{1}, GridVectors{2}, GridVectors{3});
 
    end
end

%% radiated
powerBalance.rad = calcRadiatedLossMatrix(e, h, GridVectors{1}, GridVectors{2}, GridVectors{3}, farFieldOffset);

end

function [nPorts, farFieldOffset, sigmaMasks, lumpedMasks] = parseInput(e, h, sigma, GridVectors, uMatrix, iMatrix, sMatrix, forwardPower, varInput)

% E-Field
if ~(ndims(e)==5)
    error('Invalid number of dimensions for E-Field.')
end

szFields = size(e);

if ~(szFields(end) == 3)
    error('Electric field data invalid. Needs to be of size [nPorts nx ny nz 3].');
end

if ~isfloat(e)
    error('E-Field needs to be a floating point type.');
end

% Get grid dimensions and number of ports from E-Field data
gridDimensions = szFields(2:end-1);
nPorts = szFields(1);

% H-Field
if ~(ndims(h)==5)
    error('Invalid number of dimensions for H-Field.')
end

if ~isequal(szFields, size(h))
    error('E- and H-Field do not have the same size.');
end

if ~isfloat(h)
    error('H-Field needs to be a floating point type.');
end

% sigma
if ~(ndims(sigma)==4)
    error('Invalid number of dimensions for sigma.')
end

if ~isequal([gridDimensions 3], size(sigma))
    error('Sigma size does not fit grid size.');
end

if ~isfloat(sigma)
    error('Sigma needs to be a floating point type.');
end

% grid vectors
if ~isa(GridVectors, 'cell')
    error('GridVectors needs to be a cell array of grid vectors.');
end

if ~(length(GridVectors) == 3)
    error('Three grid vectors required.');
end

if ~all(cellfun(@isvector,GridVectors));
    error('GridVectors is a cell array, but its contents are not vectors.');
end

szGrid = cellfun(@length,GridVectors);

if ~isequal(gridDimensions+1,szGrid)
    error('Grid size incompatible with field size');
end

if ~all(cellfun(@isfloat,GridVectors))
    error('Grid vectors need to be a floating point type.');
end

% U matrix
if ~ismatrix(uMatrix)
    error('Voltage matrix needs to be a matrix.');
end

if ~isfloat(uMatrix)
    error('Voltage matrix needs to be a floating point type.');
end

if ~(size(uMatrix,2) == nPorts)
    error('Voltage matrix size does not match number of ports from field data.');
end

% I matrix
if ~ismatrix(iMatrix)
    error('Current matrix needs to be a matrix.');
end

if ~isfloat(iMatrix)
    error('Current matrix needs to be a floating point type.');
end

if ~(size(iMatrix,2) == nPorts)
    error('Current matrix size does not match number of ports from field data.');
end

nLumped = size(iMatrix,1);

% S matrix
if ~ismatrix(sMatrix)
    error('Scattering matrix needs to be a matrix.');
end

if ~(size(sMatrix,1) == size(sMatrix,2))
    error('Scattering matrix needs to be square.');
end

if ~isfloat(sMatrix)
    error('Scattering matrix needs to be a floating point type.');
end

if ~(size(uMatrix,2) == nPorts)
    error('Scattering matrix size does not match number of ports from field data.');
end

% forward power
if ~(isscalar(forwardPower) && isfloat(forwardPower) && isreal(forwardPower))
    error('Forward power needs to be a floating point scalar.');
end

ip = inputParser;
addParamValue(ip, 'farFieldOffset', 4, @(x) (isnumeric(x) && isreal(x) && (fix(x) == x) && (x >= 1) ));
addParamValue(ip, 'lumpedMasks',[], @(x) (isa(x,'cell') && all(cellfun(@isvector,x)) && all(cellfun(@length,x) == nLumped) && all(cellfun(@islogical,x))));
addParamValue(ip, 'sigmaMasks',[], @(x) (isa(x,'cell') && all(cellfun(@(x) isequal(size(x), [szGrid-1 3]), x)) && all(cellfun(@islogical,x))));
parse(ip,varInput{:});

farFieldOffset = ip.Results.farFieldOffset;
sigmaMasks = ip.Results.sigmaMasks;
lumpedMasks = ip.Results.lumpedMasks;

end


