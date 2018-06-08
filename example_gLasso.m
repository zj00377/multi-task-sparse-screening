%% generate the data
n = size(X,2);
g_size = 20; %group size
Selected_feature = 5000; %The number of selected features
iteration = 10; % Iteration number
m = int16(size(X,1)*0.5); %Subsampling number

% ----------------------- set the group information --------------------- %
ind = 0:g_size:n; % the group index; the ith group include ind(i)+1:ind(i+1) features
if mod(n,g_size)~=0 
    ind(1,size(ind,2)+1)=n;
end
sg = diff(ind); % the size of each group
w = sqrt(sg)'; % the weight of each group

%----------------------- Set optional items -----------------------
opts=[];

% solver for group Lasso
opts.fName = 'glLeastR';  

% Termination 
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=1000;   % maximum number of iterations
					 % when the improvement is small, 
					 % SLEP may stop without running opts.maxIter steps
% Normalization
opts.nFlag=0;       % without normalization

% Regularization
opts.rFlag=1;       % the input parameter 'lambda' is a ratio in (0, 1]

% Group Property
opts.ind=ind;       % set the group indices
opts.q=2;           % set the value for q
%opts.sWeight=[1,1]; % set the weight for positive and negative samples
opts.gWeight=w;     % set the weight for the group, a cloumn vector
                    
opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search

%% set the regularization parameter values
% if the parameter values are the ratios of lambda/lambda_max; if use the 
% absolute value, please set opts.rFlag = 0
%
ub = 1; % upper bound of the parameter values
lb = 0.1; % lower bound of the parameter values
npar = 100; % number of parameter values
delta_lambda = (ub - lb)/(npar-1);
lambda=lb:delta_lambda:ub; % the parameter sequence
%delta_lambda = (log(ub) - log(lb))/(npar-1);
%lambda = exp(log(lb):delta_lambda:log(ub)); % the paramter sequence

%% solve the group Lasso problems along the sequence of parameter values

%tic
%[Sol, ind_zg, value] = EDPP_gLasso(X, y, lambda, opts);
%toc

%lambdaMax = max(abs(X'*y));

%for i=1:npar
%    pre(i) = size(X,2);
%    if( ind_zg(size(ind_zg,1),i)==1 )
%        pre(i) = pre(i) - (ind(1, size(ind,2))-ind(1, size(ind,2)-1));
%    end
%    pre(i) = pre(i) - nnz(ind_zg(1:size(ind_zg,1)-1, i))*g_size;
    %pre(i) = size(ind_zg, 1)-nnz(ind_zg(:,i));
%    gro(i) = nnz(Sol(:,i));
%end

%k=length(ind)-1; % the number of groups
%obj = value;

%for i=1:npar/2
%    tmp = pre(i); pre(i) = pre(npar+1-i); pre(npar+1-i)=tmp; 
%    tmp = gro(i); gro(i) = gro(npar+1-i); gro(npar+1-i)=tmp;
%    tmp = obj(i); obj(i) = obj(npar+1-i); obj(npar+1-i)=tmp;
%    tmp = lambda(i); lambda(i) = lambda(npar+1-i); lambda(npar+1-i)=tmp;
%end;

%obj(2,:) = pre;
%obj(3,:) = gro;
%obj(4,:) = lambda;
%obj = obj';

%%
t = zeros(1,size(X,2));

for i=1:iteration
    i
    k = randperm(size(X,1));
    X1 = zeros(m, size(X,2));
    y1 = zeros(m,1);
    for j=1:m/2
        X1(j,:) = X(k(j),:);  
        y1(j) = y(k(j)); 
    end
    [Sol, ind_zf] = EDPP_gLasso(X1, y1, lambda, opts);
    for j = 1:size(Sol, 1)
        t(j) = t(j)+nnz(Sol(j,:));
    end
     
    clear k x1 y1 Sol ind_zf;
end

%for i = 1:size(t, 2)
%    if( M(size(t,2)-i+1) ~= 0 )
%        s(index+1) = M(size(t,2)-i+1);
%        q(index+1) = I(size(t,2)-i+1);
%        index = index+1;
%    end
%end

[M,I] = sort(t, 'descend');
N = Selected_feature;
S = zeros(size(X,1), N);

for j=1:N
    S(:,j) = X(:,I(j));
end

X = S;

