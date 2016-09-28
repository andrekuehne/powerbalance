function pcm = calcLumpedLossMatrix(uMat,iMat)
%calcLumpedLossMatrix Calculate loss matrix from lumped element voltages and currents
%
%   pcm = calcLumpedLossMatrix(uMat, iMat)
%   calculates the power loss matrix pcm for lumped elements of an
%   n-Port device.
%
%   uMat is a nLumpedElements-by-nPorts matrix where the nth column
%   contains the lumped element voltages occurring when exciting
%   port #n.
%
%   iMat is a nLumpedElements-by-nPorts matrix where the nth column
%   contains the lumped element currents occurring when exciting
%   port #n.

pcm = (iMat'*uMat + uMat'*iMat)/4;
pcm = (pcm+pcm')/2;

end