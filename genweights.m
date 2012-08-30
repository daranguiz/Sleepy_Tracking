function weights = genweights( s, indices )
%GENWEIGHTS Generates an expected distribution given an observation.

global NUM_STATES;

% Simple
% weights = zeros(length(indices),m);
% for i=1:length(indices)
%    weights(i,indices(i)) = 1;
% end
% weights(s==0,:) = ~weights(s==0,:)+0;

% Lognormal
% distmat = abs(indices.' * ones(1,m) - ones(length(indices),1) * [1:m]);
% weights = log(distmat);
% weights = (log(s).' * ones(1,m)) - weights;
% weights = exp(-weights.^2/2/.3924^2);
% 
% % Handle the degenerate cases
% weights(distmat==0) = 0;
% weights(s==0,:) = 0;
% weights(s==0,(distmat(s==0,:)==0)) = 1;

% Gaussian
%global positions;
%weights = abs(positions(indices).' * ones(1,NUM_STATES) - ones(length(indices),1) * (1:NUM_STATES));
global sensx;
global sensy;
global statesx;
global statesy;
distx = sensx(indices)'*ones(1,NUM_STATES) - ones(length(indices),1)*statesx(:)';
disty = sensy(indices)'*ones(1,NUM_STATES) - ones(length(indices),1)*statesy(:)';
weights = sqrt(distx.^2 + disty.^2);
weights = 10./(weights.^2+1);
weights = s.' * ones(1,NUM_STATES) - weights;
weights = exp(-0.5*weights.^2);
return;

end