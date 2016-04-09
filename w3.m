function weight3 = w3(Zpt,Zpl,Iq_cartoon,patchSize,overlappingDis,sigma3)
% w3 Computes Cartoon Smoothness Weight.
%
%IN:
%	Zpt - top overlapping patch with respect to Zp
%	Zpl - left overlapping patch with respect to Zp
%	Iq_cartoon - training cartoon patch
%	patchSize - size of each patch
%	overlappingDis - overlapping distance between two adjacent patches
%   sigma3 - parameter that measures the similarity of the overlapping patches
%
%OUT:
%	weight3 - Cartoon Smoothness Weight

%{
Zpt = double(Zpt);
Zpl = double(Zpl);
Iq_cartoon = double(Iq_cartoon);

weight3 = exp(-((norm(Zpt(patchSize-overlappingDis+1:patchSize,:)-Iq_cartoon(1:overlappingDis,:))^2+...
	norm(Zpl(:,patchSize-overlappingDis+1:patchSize)-Iq_cartoon(:,1:overlappingDis))^2)/(sigma3^2)));
%}

Zpt = double(Zpt);
Zpl = double(Zpl);
Iq_cartoon = double(Iq_cartoon);

if (isempty(Zpt) && isempty(Zpl))
	weight3 = 1;
elseif (~isempty(Zpt) && isempty(Zpl))
	weight3 = exp(-(norm(Zpt(patchSize-overlappingDis+1:patchSize,:)-Iq_cartoon(1:overlappingDis,:))/...
		(sigma3^2)));
elseif (isempty(Zpt) && ~isempty(Zpl))
	weight3 = exp(-(norm(Zpl(:,patchSize-overlappingDis+1:patchSize)-Iq_cartoon(:,1:overlappingDis))/...
		(sigma3^2)));
else
	weight3 = exp(-((norm(Zpt(patchSize-overlappingDis+1:patchSize,:)-Iq_cartoon(1:overlappingDis,:))+...
		norm(Zpl(:,patchSize-overlappingDis+1:patchSize)-Iq_cartoon(:,1:overlappingDis)))/(sigma3^2)));
end
