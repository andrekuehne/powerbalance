function pcm = calcRadiatedLossMatrix(eField, hField, xLine, yLine, zLine, endOffset)
%calcMaterialLossMatrix Calculate loss matrix from power radiated into far field
%
%   Usage: pcm = calcRadiatedLossMatrix(eField, hField, sigma, xLine, yLine, zLine, endOffset)
%       eField    - multi-port electric field data in V/m on Yee grid.
%                   size(eField) should yield [nPorts nx ny nz 3]
%       hField    - multi-port magnetic field data in A/m on Yee grid.
%                   size(eField) should yield [nPorts nx ny nz 3]
%       xLine     - x coordinate vector of Yee cell vertices
%                   length(xLine) should yield nx+1
%       yLine     - y coordinate vector of Yee cell vertices
%                   length(yLine) should yield ny+1
%       zLine     - z coordinate vector of Yee cell vertices
%                   length(zLine) should yield nz+1
%       endOffset - Distance of integration box from outer
%                   boundaries in Yee cells

% get number of coils
nCoils = size(eField,1);

% define gridlines
dx=diff(xLine);
dy=diff(yLine);
dz=diff(zLine);
cxLine=xLine(1:end-1)+dx/2;
cyLine=yLine(1:end-1)+dy/2;
czLine=zLine(1:end-1)+dz/2;

% define interpolants
gEx=griddedInterpolant;
gEy=griddedInterpolant;
gEz=griddedInterpolant;
gHx=griddedInterpolant;
gHy=griddedInterpolant;
gHz=griddedInterpolant;

if nCoils == 1
    
    gEx.GridVectors={cxLine, yLine(1:end-1), zLine(1:end-1)};
    gEy.GridVectors={xLine(1:end-1), cyLine, zLine(1:end-1)};
    gEz.GridVectors={xLine(1:end-1), yLine(1:end-1), czLine};
    
    gHx.GridVectors={xLine(1:end-1), cyLine, czLine};
    gHy.GridVectors={cxLine, yLine(1:end-1), czLine};
    gHz.GridVectors={cxLine, cyLine, zLine(1:end-1)};
    
    gEx.Values=squeeze(eField(:,:,:,:,1));
    gEy.Values=squeeze(eField(:,:,:,:,2));
    gEz.Values=squeeze(eField(:,:,:,:,3));
    
    gHx.Values=squeeze(hField(:,:,:,:,1));
    gHy.Values=squeeze(hField(:,:,:,:,2));
    gHz.Values=squeeze(hField(:,:,:,:,3));
    
else
    
    gEx.GridVectors={1:nCoils, cxLine, yLine(1:end-1), zLine(1:end-1)};
    gEy.GridVectors={1:nCoils, xLine(1:end-1), cyLine, zLine(1:end-1)};
    gEz.GridVectors={1:nCoils, xLine(1:end-1), yLine(1:end-1), czLine};
    
    gHx.GridVectors={1:nCoils, xLine(1:end-1), cyLine, czLine};
    gHy.GridVectors={1:nCoils, cxLine, yLine(1:end-1), czLine};
    gHz.GridVectors={1:nCoils, cxLine, cyLine, zLine(1:end-1)};
    
    gEx.Values=eField(:,:,:,:,1);
    gEy.Values=eField(:,:,:,:,2);
    gEz.Values=eField(:,:,:,:,3);
    
    gHx.Values=hField(:,:,:,:,1);
    gHy.Values=hField(:,:,:,:,2);
    gHz.Values=hField(:,:,:,:,3);
    
end

% fill interpolants with values from simulation


if nCoils == 1
    
    % extract values
    % x directed yee faces
    xFace=xLine(end-endOffset);
    [sy, sz]=ndgrid(dy(1+endOffset:end-endOffset), dz(1+endOffset:end-endOffset));
    surfArea=sy.*sz;
    xp.surfArea=surfArea;
    xm.surfArea=-surfArea;
    
    [mx, my, mz] = ndgrid(xFace, cyLine(1+endOffset:end-endOffset), czLine(1+endOffset:end-endOffset));
    xp.ey(1,:,:,:)=gEy(mx,my,mz);
    xp.ez(1,:,:,:)=gEz(mx,my,mz);
    xp.hy(1,:,:,:)=gHy(mx,my,mz);
    xp.hz(1,:,:,:)=gHz(mx,my,mz);
    
    xFace=xLine(1+endOffset);
    [mx, my, mz] = ndgrid(xFace, cyLine(1+endOffset:end-endOffset), czLine(1+endOffset:end-endOffset));
    xm.ey(1,:,:,:)=gEy(mx,my,mz);
    xm.ez(1,:,:,:)=gEz(mx,my,mz);
    xm.hy(1,:,:,:)=gHy(mx,my,mz);
    xm.hz(1,:,:,:)=gHz(mx,my,mz);
    
    % y directed yee faces
    [sx, sz]=ndgrid(dx(1+endOffset:end-endOffset), dz(1+endOffset:end-endOffset));
    surfArea=sx.*sz;
    yp.surfArea=surfArea;
    ym.surfArea=-surfArea;
    
    yFace=yLine(end-endOffset);
    [mx, my, mz] = ndgrid(cxLine(1+endOffset:end-endOffset),yFace, czLine(1+endOffset:end-endOffset));
    yp.ex(1,:,:,:)=gEx(mx,my,mz);
    yp.ez(1,:,:,:)=gEz(mx,my,mz);
    yp.hx(1,:,:,:)=gHx(mx,my,mz);
    yp.hz(1,:,:,:)=gHz(mx,my,mz);
    
    yFace=yLine(1+endOffset);
    [mx, my, mz] = ndgrid(cxLine(1+endOffset:end-endOffset),yFace, czLine(1+endOffset:end-endOffset));
    ym.ex(1,:,:,:)=gEx(mx,my,mz);
    ym.ez(1,:,:,:)=gEz(mx,my,mz);
    ym.hx(1,:,:,:)=gHx(mx,my,mz);
    ym.hz(1,:,:,:)=gHz(mx,my,mz);
    
    % z directed yee faces
    [sx, sy]=ndgrid(dx(1+endOffset:end-endOffset), dy(1+endOffset:end-endOffset));
    surfArea=sx.*sy;
    zp.surfArea=surfArea;
    zm.surfArea=-surfArea;
    
    zFace=zLine(end-endOffset);
    [mx, my, mz] = ndgrid(cxLine(1+endOffset:end-endOffset),cyLine(1+endOffset:end-endOffset), zFace);
    zp.ex(1,:,:,:)=gEx(mx,my,mz);
    zp.ey(1,:,:,:)=gEy(mx,my,mz);
    zp.hx(1,:,:,:)=gHx(mx,my,mz);
    zp.hy(1,:,:,:)=gHy(mx,my,mz);
    
    zFace=zLine(1+endOffset);
    [mx, my, mz] = ndgrid(cxLine(1+endOffset:end-endOffset),cyLine(1+endOffset:end-endOffset), zFace);
    zm.ex(1,:,:,:)=gEx(mx,my,mz);
    zm.ey(1,:,:,:)=gEy(mx,my,mz);
    zm.hx(1,:,:,:)=gHx(mx,my,mz);
    zm.hy(1,:,:,:)=gHy(mx,my,mz);
    
else
    
    % extract values
    % x directed yee faces
    xFace=xLine(end-endOffset);
    [sy, sz]=ndgrid(dy(1+endOffset:end-endOffset), dz(1+endOffset:end-endOffset));
    surfArea=sy.*sz;
    xp.surfArea=surfArea;
    xm.surfArea=-surfArea;
    
    [mc, mx, my, mz] = ndgrid(1:nCoils, xFace, cyLine(1+endOffset:end-endOffset), czLine(1+endOffset:end-endOffset));
    xp.ey=gEy(mc,mx,my,mz);
    xp.ez=gEz(mc,mx,my,mz);
    xp.hy=gHy(mc,mx,my,mz);
    xp.hz=gHz(mc,mx,my,mz);
    
    xFace=xLine(1+endOffset);
    [mc, mx, my, mz] = ndgrid(1:nCoils, xFace, cyLine(1+endOffset:end-endOffset), czLine(1+endOffset:end-endOffset));
    xm.ey=gEy(mc,mx,my,mz);
    xm.ez=gEz(mc,mx,my,mz);
    xm.hy=gHy(mc,mx,my,mz);
    xm.hz=gHz(mc,mx,my,mz);
    
    % y directed yee faces
    [sx, sz]=ndgrid(dx(1+endOffset:end-endOffset), dz(1+endOffset:end-endOffset));
    surfArea=sx.*sz;
    yp.surfArea=surfArea;
    ym.surfArea=-surfArea;
    
    yFace=yLine(end-endOffset);
    [mc, mx, my, mz] = ndgrid(1:nCoils, cxLine(1+endOffset:end-endOffset),yFace, czLine(1+endOffset:end-endOffset));
    yp.ex=gEx(mc,mx,my,mz);
    yp.ez=gEz(mc,mx,my,mz);
    yp.hx=gHx(mc,mx,my,mz);
    yp.hz=gHz(mc,mx,my,mz);
    
    yFace=yLine(1+endOffset);
    [mc, mx, my, mz] = ndgrid(1:nCoils, cxLine(1+endOffset:end-endOffset),yFace, czLine(1+endOffset:end-endOffset));
    ym.ex=gEx(mc,mx,my,mz);
    ym.ez=gEz(mc,mx,my,mz);
    ym.hx=gHx(mc,mx,my,mz);
    ym.hz=gHz(mc,mx,my,mz);
    
    % z directed yee faces
    [sx, sy]=ndgrid(dx(1+endOffset:end-endOffset), dy(1+endOffset:end-endOffset));
    surfArea=sx.*sy;
    zp.surfArea=surfArea;
    zm.surfArea=-surfArea;
    
    zFace=zLine(end-endOffset);
    [mc, mx, my, mz] = ndgrid(1:nCoils, cxLine(1+endOffset:end-endOffset),cyLine(1+endOffset:end-endOffset), zFace);
    zp.ex=gEx(mc,mx,my,mz);
    zp.ey=gEy(mc,mx,my,mz);
    zp.hx=gHx(mc,mx,my,mz);
    zp.hy=gHy(mc,mx,my,mz);
    
    zFace=zLine(1+endOffset);
    [mc, mx, my, mz] = ndgrid(1:nCoils, cxLine(1+endOffset:end-endOffset),cyLine(1+endOffset:end-endOffset), zFace);
    zm.ex=gEx(mc,mx,my,mz);
    zm.ey=gEy(mc,mx,my,mz);
    zm.hx=gHx(mc,mx,my,mz);
    zm.hy=gHy(mc,mx,my,mz);
end
% construct matrices
% x directed yee faces
xp.surfArea = permute(xp.surfArea(:),[2 1]);
xp.ey = reshape(xp.ey,nCoils,[]);
xp.ez = reshape(xp.ez,nCoils,[]);
xp.hy = reshape(xp.hy,nCoils,[]);
xp.hz = reshape(xp.hz,nCoils,[]);
mxp = bsxfun(@times,xp.ey,xp.surfArea)*xp.hz'-bsxfun(@times,xp.ez,xp.surfArea)*xp.hy';

xm.surfArea = permute(xm.surfArea(:),[2 1]);
xm.ey = reshape(xm.ey,nCoils,[]);
xm.ez = reshape(xm.ez,nCoils,[]);
xm.hy = reshape(xm.hy,nCoils,[]);
xm.hz = reshape(xm.hz,nCoils,[]);
mxm = bsxfun(@times,xm.ey,xm.surfArea)*xm.hz'-bsxfun(@times,xm.ez,xm.surfArea)*xm.hy';

% y directed yee faces
yp.surfArea = permute(yp.surfArea(:),[2 1]);
yp.ex = reshape(yp.ex,nCoils,[]);
yp.ez = reshape(yp.ez,nCoils,[]);
yp.hx = reshape(yp.hx,nCoils,[]);
yp.hz = reshape(yp.hz,nCoils,[]);
myp = bsxfun(@times,yp.ez,yp.surfArea)*yp.hx'-bsxfun(@times,yp.ex,yp.surfArea)*yp.hz';

ym.surfArea = permute(ym.surfArea(:),[2 1]);
ym.ex = reshape(ym.ex,nCoils,[]);
ym.ez = reshape(ym.ez,nCoils,[]);
ym.hx = reshape(ym.hx,nCoils,[]);
ym.hz = reshape(ym.hz,nCoils,[]);
mym = bsxfun(@times,ym.ez,ym.surfArea)*ym.hx'-bsxfun(@times,ym.ex,ym.surfArea)*ym.hz';

% z directed yee faces
zp.surfArea = permute(zp.surfArea(:),[2 1]);
zp.ex = reshape(zp.ex,nCoils,[]);
zp.ey = reshape(zp.ey,nCoils,[]);
zp.hx = reshape(zp.hx,nCoils,[]);
zp.hy = reshape(zp.hy,nCoils,[]);
mzp = bsxfun(@times,zp.ex,zp.surfArea)*zp.hy'-bsxfun(@times,zp.ey,zp.surfArea)*zp.hx';

zm.surfArea = permute(zm.surfArea(:),[2 1]);
zm.ex = reshape(zm.ex,nCoils,[]);
zm.ey = reshape(zm.ey,nCoils,[]);
zm.hx = reshape(zm.hx,nCoils,[]);
zm.hy = reshape(zm.hy,nCoils,[]);
mzm = bsxfun(@times,zm.ex,zm.surfArea)*zm.hy'-bsxfun(@times,zm.ey,zm.surfArea)*zm.hx';

% sum values
pcm = 0.5*(mxp+mxm+myp+mym+mzp+mzm).';
pcm = (pcm+pcm')/2;



end
