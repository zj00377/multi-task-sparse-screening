

%% set up the solver from SLEP, can leave as default
opts=[];

% termination criterion
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=1000;   % maximum number of iterations
					 % when the improvement is small, 
					 % SLEP may stop without running opts.maxIter steps

% normalization
opts.nFlag=0;       % without normalization

% regularization
opts.rFlag=1;       % the input parameter 'lambda' is a ratio in (0, 1]

opts.fName = 'LeastR';  % compute a sequence of lasso problems

opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search

%% set the regularization paramter values
% if the parameter values are the ratios of lambda/lambda_max; if use the 
% absolute value, please set opts.rFlag = 0
%
ub = 1.0; % upper bound of the parameter values
lb = 0.1; % lower bound of the parameter values
npar = 100; % number of parameter values
delta_lambda = (ub - lb)/(npar-1);
lambda=lb:delta_lambda:ub; % the parameter sequence
%delta_lambda = (log(ub) - log(lb))/(npar-1);
%lambda = exp(log(lb):delta_lambda:log(ub)); % the paramter sequence

%% solve the Lasso problems along the sequence of parameter values

tic
[Sol, ind_zf] = EDPP_Lasso(X, y, lambda, opts);
toc

lambdaMax = max(abs(X'*y));

for i=1:npar
    pre(i) = size(X,2)-nnz(ind_zf(:,i)); 
    gro(i) = nnz(Sol(:,i));
end

for i=1:npar
    obj(i) = 0.5*sum((X*Sol(:,i)-y).^2)+lambda(i)*lambdaMax*sum(abs(Sol(:,i))); 
end;

for i=1:npar/2
    tmp = pre(i); pre(i) = pre(npar+1-i); pre(npar+1-i)=tmp; 
    tmp = gro(i); gro(i) = gro(npar+1-i); gro(npar+1-i)=tmp;
    tmp = obj(i); obj(i) = obj(npar+1-i); obj(npar+1-i)=tmp;
    tmp = lambda(i); lambda(i) = lambda(npar+1-i); lambda(npar+1-i)=tmp;
end;

obj(2,:) = pre;
obj(3,:) = gro;
obj(4,:) = lambda;
obj = obj';



