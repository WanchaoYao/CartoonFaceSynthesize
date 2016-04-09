function weight1 = w1(Xp,Iq_photo,sigma1)
% w1 Computes Photo-Photo Weight.
%
%IN:
%	Xp - input photo patch
%	Iq_photo - training photo patch
%   sigma1 - parameter that adjusts the range(i.e., intensity) similarity
%
%OUT:
%	weight1 - Photo-Photo Weight

Xp = double(Xp);
Iq_photo = double(Iq_photo);

Xp_normalized = Xp-mean(Xp(:));
Iq_photo_normalized = Iq_photo-mean(Iq_photo(:));

weight1 = exp(-((norm(Xp_normalized-Iq_photo_normalized)^2)/(sigma1^2)));