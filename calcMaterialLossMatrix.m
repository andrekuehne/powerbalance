function pcm = calcMaterialLossMatrix(eField, sigma, xLine, yLine, zLine)
%calcMaterialLossMatrix Calculate loss matrix from electric fields in lossy materials
%
%   Usage: pcm = calcMaterialLossMatrix(eField, sigma, xLine, yLine, zLine)
%       eField  - multi-port electric field data in V/m on Yee grid.
%                 size(eField) should yield [nPorts nx ny nz 3]
%       sigma   - Electric conductivity in S/m
%                 size(sigma) should yield [nx ny nz 3]
%       xLine   - x coordinates of Yee cell vertices
%                 length(xLine) should yield nx+1
%       yLine   - y coordinates of Yee cell vertices
%                 length(yLine) should yield ny+1
%       zLine   - z coordinates of Yee cell vertices
%                 length(zLine) should yield nz+1

% generate grid variables
dx=diff(xLine);
dy=diff(yLine);
dz=diff(zLine);
cxLine=xLine(1:end-1)+dx/2;
cyLine=yLine(1:end-1)+dy/2;
czLine=zLine(1:end-1)+dz/2;
dcx=diff(cxLine);
dcy=diff(cyLine);
dcz=diff(czLine);

% generate grids for integration deltas
[mxx, myx, mzx] = ndgrid(dx, dcy ,dcz);
vx=mxx.*myx.*mzx;
[mxy, myy, mzy] = ndgrid(dcx, dy ,dcz);
vy=mxy.*myy.*mzy;
[mxz, myz, mzz] = ndgrid(dcx, dcy ,dz);
vz=mxz.*myz.*mzz;

% pad grid with zeros for unusable boundary regions
vx = padarray(vx,[0 1 1],0,'pre');
vy = padarray(vy,[1 0 1],0,'pre');
vz = padarray(vz,[1 1 0],0,'pre');

% concatenate to full 3D vector grid
v = cat(4,vx,vy,vz);
clear vx vy vz;

% find points with nonzero conductivity
validEdges = sigma>0;

% only use cell edges with conductivity
eValid = eField(:,validEdges);
sValid = sigma(validEdges)';
vValid = v(validEdges)';

% calculate power correlation matrix
pcm = 0.5*(bsxfun(@times,eValid,sValid.*vValid)*eValid').';
pcm=(pcm+pcm')/2; %enforce hermiticity



end