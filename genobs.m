function s = genobs( b, indices )
%GENOBS Generates observations for the given indices of positions, and
%using b as the actual object location.

% Simple
%s = (indices==b)+0;

% Lognormal
% s = zeros(1,length(indices));
% temp = indices(indices~=b);
% s(indices~=b) = exp( log(abs(temp-b)) + 0.3924*randn(size(temp)) );

% Gaussian
%global positions;
global sensx;
global sensy;
global statesx;
global statesy;
dist = sqrt((sensx(indices) - statesx(b)).^2 + (sensy(indices) - statesy(b)).^2);
%s = 10./( (positions(indices)-b).^2 + 1 ) + randn(1,length(indices));
s = 10./( dist.^2 + 1 ) + 0.5*randn(1,length(indices));

return;

end