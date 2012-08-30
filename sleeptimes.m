function sleepvals=sleeptimes( dist, indices, initflag )

global AWAKE_COST;
c = AWAKE_COST;
global P;

global V;
global Vinv;
global lambdas;
global rotational;
global redmat;

global T;

global J;
global uvals;

global NUM_STATES;
m = NUM_STATES;
global NUM_SENSORS;
n = NUM_SENSORS;

% The first time this function is called with a new P matrix, it should be
% called with initflag equal to true so that P can be diagonalized and the
% policy computed.
if (initflag ~= 0)
    redind = (1:m).';

    % We diagonalize P so that P = V * diag(lambdas) * Vinv.  We order the
    % eigenvalues from largest to smallest.  Also, if period > 1, we also
    % group eigenvalues whose periodth powers are equal.
    [V,temp] = eig( full(P) );
    lambdas = temp(logical(eye(size(temp))));
    
    temp = lambdas;
    notused = (1:m).';
    order = zeros(m,1);
    for i=0:m-1
        if (rem(i,1)==0)
            [~,pos] = max(abs(temp(notused)));
        else
            [~,pos] = min(abs(temp(notused)-last));
        end
        pos = notused(pos);
        order(i+1) = pos;
        last = temp(pos);
        notused = notused(notused~=pos);
    end
    
    lambdas = lambdas(order);
    V = V(:,order);
    Vinv = V\eye(size(P));
    

    % When computing our upper and lower bounds, we need to know whether the
    % eigenvalues we are using are nonnegative real numbers.  If not, then
    % these eigenvalues will rotate as we increase u so we call them
    % rotational.
    temp = lambdas(redind);
    rotational = find(abs(imag(temp)) > 1e-4*real(temp));

    % Define a reduction matrix that is designed to consolidate values that
    % correspond to eigenvalues of the same magnitude.
    redmat = spalloc( m, ceil(m), m );
    for i=1:ceil(m)
        redmat((i:min(i,m)),i) = 1;
    end
    
    % Now perform policy iteration to define the QMDP sleeping policy.
    
    % The initial policy is to always sleep forever.
    if ( any(size(J)~=[m,n]) || any(size(uvals)~=[m,n]) )
        J = (eye(m)-P) \ T;
        uvals = Inf * ones(m,n);
    end
    tcosts = zeros(m,n);     % tracking costs
    ecosts = zeros(m,n);     % energy costs

    % Define some variables used in the iterations.
    newuvals = zeros(size(uvals));
    A = zeros(m);
    tx = zeros(m,1);
    ex = zeros(m,1);
    invmat = inv(eye(m)-P);

    flag = 0;
    iter = 0;
    while (~flag)
        iter = iter+1;
        prevJsum = sum(sum(J));

        % Find a new matrix of sleep times.
        for state=1:m
            vec = [zeros(1,state-1) 1 zeros(1,m-state)];
            newuvals(state,:) = dowork( vec, (1:n), 0 );
        end

        % Use the new sleep times to update the J matrix.  This involves
        % solving a set of linear equations for each sensor.
        for sensor=1:n
            for state=1:m
                u = newuvals(state,sensor);
                if (u == Inf)
                    A(state,:) = 0;
                    tx(state) = invmat(state,:) * T(:,sensor);
                    ex(state) = 0;
                else
                    temp = real((V(state,:) .* (lambdas.^u).') * Vinv);
                    A(state,:) = temp*P;
		            tx(state) = (invmat(state,:) - temp * invmat) * T(:,sensor);
                    ex(state) = c*sum(temp*P);
                end
            end
            A = eye(size(A)) - A;
            tcosts(:,sensor) = A \ tx;
            ecosts(:,sensor) = A \ ex;
        end
        J = tcosts + ecosts;

        % Terminate if the matrix of sleep times did not change.
        if ( ~any( newuvals ~= uvals ) )
            flag = 1;
        end
        if ( (sum(sum(J)) > prevJsum) )
            flag = 1;
        end
        if (iter == 1)
            flag = 0;
        end
        max(newuvals(newuvals~=Inf))
        uvals = newuvals;
    end
%    normfact = sum(dist*P*invmat);
%    x2 = [x2 ; sum(dist*ecosts)/c/normfact];
%    y2 = [y2 ; sum(dist*tcosts)/normfact];
end

% Actually compute the sleep times.
sleepvals = dowork( dist, indices, 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=dowork( dist, indices, trainingdone )

global AWAKE_COST;
c = AWAKE_COST;
global P;

global V;
global Vinv;
global lambdas;
global rotational;
global redmat;

global T;

global J;
global uvals;

global NUM_STATES;
m = NUM_STATES;

% If the input distribution is a point mass distribution, the sleep times
% can be computed via table lookup.
if (trainingdone && (sum(dist > 0) == 1))
    y = uvals( dist > 0, indices );
    return;
end

numind = length(indices);

redind = (1:m).';

y = Inf * ones(size(indices));

% Convert the distribution to its diagonalized form
dist = dist*V;

% The ainit matrix is a numind by m matrix that when multiplied
% by lambdas.^u yields the function we are trying to minimize evaluated
% for each index at a sleep time of u.
ainit = (Vinv*P*(c + J(:,indices))).';
ainit = ainit - (diag(1./(1-lambdas)) * Vinv * T(:,indices)).';
ainit = (ones(numind,1)*dist) .* ainit;

L = sparse(diag(lambdas(redind)));

mincosts = zeros(numind,1);

% If there are multiple maximal eigenvalues, we successively consider the
% cases when u mod period = iter for iter=0 to period-1.
for iter=0:0
    u = iter;
    avals = ainit * redmat;
	ndlist = (1:numind).';  % The list of "not done" indices
	while ( ~isempty(ndlist) )
        % Compute the function for a sleep time of u.
        acosts = real( sum(avals,2) );
        
        % Compute a lower bound on our function for sleep times >= u.
        lbtemp = avals;
        lbtemp(:,rotational) = -abs(lbtemp(:,rotational));
        lbtemp = real( lbtemp );
               
        for j=1:size(lbtemp,2)-1
            lbtemp(:,j+1) = lbtemp(:,j+1) + max(lbtemp(:,j),0);
        end
        lbtemp = min(lbtemp,0);

        lb = sum(lbtemp,2);
        
        % Update y if we have found any new minima.        
        temp = find(acosts < mincosts(ndlist));
        mincosts(ndlist(temp)) = acosts(temp);
        temp = ndlist(temp);
        y(temp) = u;
        
        % See if we have found the global minima for any of our indices.
        temp = (lb >= mincosts(ndlist));
        ndlist = ndlist(~temp);
        avals = avals(~temp,:);

        % Prepare for the next iteration.
        u = u+1;
        avals = avals * L;
	end
	ainit = ainit .* (ones(numind,1)*(lambdas.'));
end
