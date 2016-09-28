function pcm = calcForwardPowerMatrix(nPorts, forwardPower)
%calcForwardPowerMatrix Calculate forward power matrix from forward power and number of device ports
%   nPorts - number of ports of device
%   forwardPower - Forward power incident into ports

pcm = eye(nPorts) * forwardPower;

end
