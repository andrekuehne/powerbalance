%% load exported data
% SParameter import via S-Parameter-toolbox from the Matlab file exchange:
% https://de.mathworks.com/matlabcentral/fileexchange/6080-s-parameter-toolbox--+-z--y--h--g--abcd--t-

load '2ch_loops.exported/e_field.mat';
load '2ch_loops.exported/h_field.mat';
load '2ch_loops.exported/current_density.mat';
load '2ch_loops.exported/u_matrix.mat';
load '2ch_loops.exported/i_matrix.mat';
load '2ch_loops.exported/meshdata.mat';
[f,sMat] = SXPParse('2ch_loops.exported/s_matrix.s2p');

%% crop last points in field data
e = e(:,1:end-1,1:end-1,1:end-1,:);
h = h(:,1:end-1,1:end-1,1:end-1,:);
c = c(:,1:end-1,1:end-1,1:end-1,:);

%% calculate conductivity
sigma = squeeze(real(c(1,:,:,:,:)./e(1,:,:,:,:)));
sigma(~isfinite(sigma)) = 0;

%% generate binary masks for differentiation of different materials
maskFR4 = sigma>0 & sigma < 0.001;
maskPhantom = sigma>0 & sigma >= 0.001;

%% generate binary masks for lumped component differentiation
maskCM = logical([1 1 0 0 0 0 0 0 0]); % matching caps
maskCT = logical([0 0 1 1 1 1 1 1 0]); % fixed & tuning caps
maskCD = logical([0 0 0 0 0 0 0 0 1]); % decoupling cap

%% calculate power balance
pBal = calcPowerBalance(e,h,sigma,{xLine,yLine,zLine},uMat,iMat,sMat,0.02,'sigmaMasks',{maskFR4,maskPhantom},'lumpedMasks',{maskCM,maskCT,maskCD},'farfieldOffset',2);

%% example plot
nPoints = 361;
phs = linspace(0,180,nPoints);
v = [ones(1,nPoints); exp(1j*phs/180*pi)]./sqrt(2);

pRef = zeros(nPoints,1);
pRad = zeros(nPoints,1);
pPhantom = zeros(nPoints,1);
pSubstrate = zeros(nPoints,1);
pMatch = zeros(nPoints,1);
pFixed = zeros(nPoints,1);
pDecouple = zeros(nPoints,1);

for k=1:nPoints
    pRef(k) = real(v(:,k)'*pBal.cpl*v(:,k));
    pRad(k) = real(v(:,k)'*pBal.rad*v(:,k));
    pPhantom(k) = real(v(:,k)'*pBal.matMasked{2}*v(:,k));
    pSubstrate(k) = real(v(:,k)'*pBal.matMasked{1}*v(:,k));
    pMatch(k) = real(v(:,k)'*pBal.lmpMasked{1}*v(:,k));
    pFixed(k) = real(v(:,k)'*pBal.lmpMasked{2}*v(:,k));
    pDecouple(k) = real(v(:,k)'*pBal.lmpMasked{3}*v(:,k));
end

figure('color','w');
area(phs,[pPhantom,pSubstrate,pMatch,pFixed,pDecouple,pRad,pRef]./0.02*100);
axis([0 180 0 100]);
xlabel('Relative channel phase [°]');
ylabel('Fraction of forward power [%]');
legend('Phantom','Substrate','C_{Match}','C_{Tune}','C_{Decouple}','Radiated','Reflected','Location','SouthWest');
title('Power balance of a 2-channel loop array under phase-only shimming');