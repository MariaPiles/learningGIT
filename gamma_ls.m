%
% gamma_ls.m
% [out,W,bestmu,M,nMSE]=gamma_ls(x,y,xtst,ytst,dibuja)
%
% This function adapts a gamma filter using least squares (LS)
% 
% Inputs:		
%		x:      input sequence
%		y:      (desired) output sequence
%       xtst:   test input sequence
%       ytst:   test output sequence
%       dibuja: flag for plotting the optimization
%
% Outputs:
%		W:      weight vector
%		bestmu: best mu parameter between 0 and 2
%		M:      memory depth = K/bestmu; 
% 
%
% (c)   Gus, 2002, e-mail:  gustavo.camps@uv.es
%       Mofified: JL        jlrojo@tsc.uc3m.es

function [W,mubest,M]=gamma_ls(x,y,xtst,ytst,dibuja)

% Iniciales
if nargin==2, dibuja=0; end
y=y(:); x=x(:);
N = length(y);
mus    = 0.01:0.01:.99;
kas    = 1:30;
kbest  = 10;
mubest = .25;

% the iterative hyperparam tuning trick ...
for i=1:5  % enough!
    
    % Barrido en mu
    errmu=zeros(size(mus));
    for m=1:length(mus)
        W = gammasolve(x,y,mus(m),kbest); 
        otst=salidagamma(xtst,W,mus(m));
        e=otst(length(W)+1:end)-ytst(length(W)+1:end);
        errmu(m) = sum(e.^2)/length(e);
    end
    [~,m] = min(errmu);
    mubest=mus(m);
    
    % Barrido en k
    errk=zeros(size(kas));
    for m=1:length(kas)
        W = gammasolve(x,y,mubest,kas(m)); 
        otst=salidagamma(xtst,W,mubest);
        e=otst(length(W)+1:end)-ytst(length(W)+1:end);
        errk(m) = sum(e.^2)/length(e);
    end
    [~,m] = min(errk);
    kbest=kas(m);
    
    if dibuja
        figure(1)
        subplot(211), plot(mus,errmu); axis tight; xlabel('\mu');ylabel('MSE'),hold on
        subplot(212), plot(kas,errk), axis tight, xlabel('K');ylabel('MSE')
        drawnow,
    end
end


% The best output is:
[W,~,~] = gammasolve(x,y,mubest,kbest);
M = kbest/mubest;


% -------------------------------------------------
% ---------------- Auxiliars ---------------------- 
% -------------------------------------------------

function ypred=salidagamma(x,W,mu)

% Construct X matrix and w vector  
p=length(W);
k0 = p+1;
X=[];
bgamma=mu;
agamma=[1 -(1-mu)];
xold=x';
for i=1:p
    X=[X;xold];
    xold=filter(bgamma,agamma,xold);
end
X=X(:,k0:end);

ypred=X'*W;
ypred=[zeros(p,1);ypred];
ypred=ypred(:);


function [W,err,ypred] = gammasolve(x,y,mu,k)

% Construct X matrix and pseudoinversion  
X=[];
bgamma=mu;
agamma=[1 -(1-mu)];
xold=x';
for i=1:k
    X=[X;xold];
    xold=filter(bgamma,agamma,xold);
end
X=X(:,k+1:end);

yr=y(k+1:end);
W=pinv(X')*yr;
ypred=X'*W;
err=sum( (ypred-yr).^2)/length(yr);
