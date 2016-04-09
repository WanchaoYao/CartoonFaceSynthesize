function weight4 = w4(Cp,Cq,sigma4)
% w4 Computes Spatial Weight.
%
%IN:
%	Cp - center coordinate of patch p
%	Cq - center coordinate of patch q
%   sigma4 - parameter that adjusts the spatial similarity
%
%OUT:
%	weight4 - Spatial Weight

Cp = double(Cp);
Cq = double(Cq);

weight4 = exp(-((norm(Cp-Cq)^2)/(sigma4^2)));