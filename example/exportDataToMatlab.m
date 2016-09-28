%% Load CST library
try 
    unloadlibrary('CSTResultReader_AMD64');
catch
end

DLL_Path = '..\';
HeaderPath = 'C:\Program Files (x86)\CST STUDIO SUITE 2016\ResultReaderDLL\';
libname = 'CSTResultReader';
if (strcmp(computer, 'PCWIN64'))
    libname = 'CSTResultReader_AMD64';
    DLL_Path = ['C:\Program Files (x86)\CST STUDIO SUITE 2016\AMD64\'];
end
CSTResultReader = strcat( DLL_Path, libname, '.dll' );
CSTResultReaderH = strcat( HeaderPath, 'CSTResultReaderInterf.h' );
loadlibrary(CSTResultReader,CSTResultReaderH);
sHandle.m_pProj = 0;

%% initialize project
fullPathName = '2ch_loops.cst';
[ret, mwsProjName, sHandle] = calllib(libname, 'CST_OpenProject', fullPathName, sHandle);  % Initialize the handle

%% read mesh data
nxyz = zeros(1,3,'int32');                                                              % Initialize variable for mesh dimensions
[ret, sHandle, nxyz] = calllib(libname, 'CST_GetHexMeshInfo', sHandle, nxyz); % Read out hex mesh dimensions
nx = nxyz(1); ny = nxyz(2); nz = nxyz(3);                                               % Store mesh dimensions separately
% Constants for reparametrization/super-index applied later:
% n = 1 + (i-1) + (j-1)*nx + (k-1)*nx*ny
nmax = nx*ny*nz;                % maximum value of super-index n
nCells = (nx-1)*(ny-1)*(nz-1);  % number of cells

xyzLines = zeros(1,sum(nxyz),'double');                                                         % Initialize variable for mesh lines
[ret, sHandle, xyzLines] = calllib(libname, 'CST_GetHexMesh', sHandle, xyzLines);     % Read out hex mesh lines
xLine = xyzLines(1:nx); yLine = xyzLines(nx+1:nx+ny); zLine = xyzLines(nx+ny+1:nx+ny+nz);    % Store mesh line positions separately

nVoxels=numel(xLine)*numel(yLine)*numel(zLine);
save('2ch_loops.exported/meshdata.mat','xLine','yLine','zLine');
%% get E-Fields
freq = '297.2';
portPrefix = 'SS_TxCh';
nPorts = 2;

iResultNumber = 0;
e = zeros([nPorts 3*nVoxels]);
for fieldCnt=1:nPorts
    

    sTree3DResName = sprintf('2D/3D Results\\E-Field\\e-field (f=%s) [%s%d]', freq, portPrefix, fieldCnt);

    resSize=0;
    % Ask for the size of the expected data vector.
    [ret, sHandle, sTree3DResName, resSize] = calllib(libname, 'CST_Get3DHexResultSize', sHandle, sTree3DResName, iResultNumber, resSize);
    assert(~ret)
    
    % Define an array of proper size.
    % NOTE: If the array size does not have the proper size, matlab may crash!
    fTmp=zeros(resSize,1,'single');
    
    % Get the e-field values.
    [ret, sHandle, sTree3DResName, fTmp] = calllib(libname, 'CST_Get3DHexResult', sHandle, sTree3DResName, iResultNumber, fTmp);
    assert(~ret)
    
    e(fieldCnt,:) = complex(fTmp(1:2:end),fTmp(2:2:end));
    
end
e = (reshape(e, [nPorts numel(xLine) numel(yLine) numel(zLine) 3]));
save('2ch_loops.exported/e_field.mat','e');
clear e;

%% get H-Fields
freq = '297.2';
portPrefix = 'SS_TxCh';
nPorts = 2;

iResultNumber = 0;
h = zeros([nPorts 3*nVoxels]);
for fieldCnt=1:nPorts
    

    sTree3DResName = sprintf('2D/3D Results\\H-Field\\h-field (f=%s) [%s%d]', freq, portPrefix, fieldCnt);

    resSize=0;
    % Ask for the size of the expected data vector.
    [ret, sHandle, sTree3DResName, resSize] = calllib(libname, 'CST_Get3DHexResultSize', sHandle, sTree3DResName, iResultNumber, resSize);
    assert(~ret)
    
    % Define an array of proper size.
    % NOTE: If the array size does not have the proper size, matlab may crash!
    fTmp=zeros(resSize,1,'single');
    
    % Get the e-field values.
    [ret, sHandle, sTree3DResName, fTmp] = calllib(libname, 'CST_Get3DHexResult', sHandle, sTree3DResName, iResultNumber, fTmp);
    assert(~ret)
    
    h(fieldCnt,:) = complex(fTmp(1:2:end),fTmp(2:2:end));
    
end
h = (reshape(h, [nPorts numel(xLine) numel(yLine) numel(zLine) 3]));
save('2ch_loops.exported/h_field.mat','h');
clear h;

%% get current density
freq = '297.2';
portPrefix = 'SS_TxCh';
nPorts = 2;

iResultNumber = 0;
c = zeros([nPorts 3*nVoxels]);
for fieldCnt=1:nPorts
    

    sTree3DResName = sprintf('2D/3D Results\\Current Density\\current (f=%s) [%s%d]', freq, portPrefix, fieldCnt);

    resSize=0;
    % Ask for the size of the expected data vector.
    [ret, sHandle, sTree3DResName, resSize] = calllib(libname, 'CST_Get3DHexResultSize', sHandle, sTree3DResName, iResultNumber, resSize);
    assert(~ret)
    
    % Define an array of proper size.
    % NOTE: If the array size does not have the proper size, matlab may crash!
    fTmp=zeros(resSize,1,'single');
    
    % Get the e-field values.
    [ret, sHandle, sTree3DResName, fTmp] = calllib(libname, 'CST_Get3DHexResult', sHandle, sTree3DResName, iResultNumber, fTmp);
    assert(~ret)
    
    c(fieldCnt,:) = complex(fTmp(1:2:end),fTmp(2:2:end));
    
end
c = (reshape(c, [nPorts numel(xLine) numel(yLine) numel(zLine) 3]));
save('2ch_loops.exported/current_density.mat','c');
clear c;

%% get probe voltages
sTree1DName = '1D Results\DS_Results\Lumped_Elements\Voltages';

% Specify the desired result number. In most cases there is only one result.
% (as assumed here)
iResultNumber = 0;

% Size of expected 1D-Data
nSigSize = 0;

% Ask for the size of the expected data vector.
[ret, sHandle, sTree1DName, nSigSize] = ...
    calllib(libname, 'CST_Get1DResultSize', sHandle, sTree1DName, iResultNumber, nSigSize);

xVal=zeros(nSigSize,1);

% Get the x axis (Tx channel number) values.
[ret, sHandle, sTree1DName, xVal] = ...
    calllib(libname, 'CST_Get1DRealDataAbszissa', sHandle, sTree1DName, iResultNumber, xVal);

% Get the complex voltage values.
yVal=zeros(nSigSize*2,1);
[ret, sHandle, sTree1DName, yVal] = ...
    calllib(libname, 'CST_Get1D_2Comp_DataOrdinate', sHandle, sTree1DName, iResultNumber, yVal);
voltages = complex(yVal(1:2:end),yVal(2:2:end));

nTxCh = max(xVal);
nLumped = numel(xVal)/nTxCh;
uMat = reshape(voltages,nLumped,nTxCh);
save('2ch_loops.exported/u_matrix.mat','uMat');

%% get probe currents
sTree1DName = '1D Results\DS_Results\Lumped_Elements\Currents';

% Specify the desired result number. In most cases there is only one result.
% (as assumed here)
iResultNumber = 0;

% Size of expected 1D-Data
nSigSize = 0;

% Ask for the size of the expected data vector.
[ret, sHandle, sTree1DName, nSigSize] = ...
    calllib(libname, 'CST_Get1DResultSize', sHandle, sTree1DName, iResultNumber, nSigSize);

xVal=zeros(nSigSize,1);

% Get the x axis (Tx channel number) values.
[ret, sHandle, sTree1DName, xVal] = ...
    calllib(libname, 'CST_Get1DRealDataAbszissa', sHandle, sTree1DName, iResultNumber, xVal);

% Get the complex current values.
yVal=zeros(nSigSize*2,1);
[ret, sHandle, sTree1DName, yVal] = ...
    calllib(libname, 'CST_Get1D_2Comp_DataOrdinate', sHandle, sTree1DName, iResultNumber, yVal);
currents = complex(yVal(1:2:end),yVal(2:2:end));

nTxCh = max(xVal);
nLumped = numel(xVal)/nTxCh;
iMat = reshape(currents,nLumped,nTxCh);
save('2ch_loops.exported/i_matrix.mat','iMat');