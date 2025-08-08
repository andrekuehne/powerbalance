 bnfunction pcm = calcCouplingLossMatrix(sMat, forwardPower)
%calcCouplingLossMatrix Calculate coupling loss matrix from scattering matrix and forward power
%   sMat - scattering matrix
%   forwardPower - Forward power incident into ports

pcm = (sMat'*sMat) * forwardPower;

end