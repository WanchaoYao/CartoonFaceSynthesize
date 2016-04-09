function weight2 = w2(Dp,Dq_cartoon,sigma2)
% w2 Computes Photo-Cartoon Weight.
%
%IN:
%	Dp - PHOG of input photo patch
%	Dq_cartoon - PHOG of training cartoon patch
%   sigma2 - parameter that adjusts the descriptor similarity
%
%OUT:
%	weight2 - Photo-Cartoon Weight

Dp = double(Dp);
Dq_cartoon = double(Dq_cartoon);

weight2 = exp(-((norm(Dp-Dq_cartoon)^2)/(sigma2^2)));