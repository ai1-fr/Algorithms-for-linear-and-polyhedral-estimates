function [dta,cntr] = set_ctl(n,K)
% [DATA,CONTROL] = SET_CTL
% Set DATA and CONTROL structures for CTL experiments with random A matrix
%   By A. Juditsky 03/31/2024
%

ni=nargin;
if ni<2, K=8; end
if ~ni, n=128; end

seed=144;
cntr.n=n; % signal dimension
cntr.tolbase=1; % define basic numeric tolerance parameter
cntr.S=10;  % number of pieces in the bundle, parameter tau in the paper
cntr.I=10;  % number of components in a piece, parameter rho in the paper
cntr.tol=1.1; % required relative gap
cntr.bmode='P'; % bundle mode setting: 'P'=plain, 'A'=advanced
cntr.pi=0; % pi=1 for Euclidean algorithm, pi=0 for Mirror Descent
cntr.sol='M'; % IPM solver to be used by CVX, 'M'=mosek, 'S'=SDPT3
cntr.print=10000; % control output convergence information
                  % when PRINT>=10000, output no information
cntr.R=100; % norm bound
cntr.itrmax=3000; % max #oracle calls


% if nargin, dta.ID=id; else dta.ID=1; end % experiment identifier; 
%                                            % save results in ['data',num2str{data.ID),'.mat'] file
dta.n=cntr.n;
dta.stype=1; % signal set type: 1=intersection of co-axial elliptic cylinders
               % 2=smooth signal
dta.box=0; % when SigType=I, setting dta.box=1 will result in X which is intesection of `Sobolev parallelepipeds', 
           % when dta.box=0, X is intersection of Sobolev ellipsoids
dta.K=K; % number of intersected elliptic cylinders when X is the intersection of cylinders
          % or number of discretizations for smooth signal
dta.dec=0.5; % exponent of the decay of the axes of elliptic cylinders
dta=create_sT(dta);
dta.p=inf; % ellitope T is the unit ball of the p-norm of blocks of size n/K, 1<=p<=inf 
if dta.p<inf 
    dta.q=dta.p/(dta.p-1);
else
    dta.q=1;
end
dta.m=dta.n; % square matrix A
% dta.I=cntr.I;   
% define matrix parameteres A and B here 
dta.Atype='G'; % A has iid random entries of Gaussian ('G') or Rademacher ('R') type
               % or it can be randomly rotated matrix ('O') with singular values decaying with exponent dta.edr 
dta.edr=1; % decay exponent of the singular values of A of type 'O''
dta.Btype='S'; % matrix B defining recovery metric:'D' for recovery of derivative, 'P' primitive, 
               % 'S' for identity matrix (signal recovery) 
dta.sigma=0.1; % noise sigma
dta.eps=0.05; % desired reliability of the estimate
dta=create_AB(dta,seed); % create matrices A and B

cntr.K=dta.K; % same as dta.K 
end

function dta=create_sT(dta)
FCTR=0.2;
K=dta.K;
%
dn=ceil(dta.n/dta.K);
stype=dta.stype;
if stype==1 % ellitopic X
    xMax=zeros(dta.n,1);
    if dta.box
        dta.sT=cell(dta.K,1);
        dta.dgsT=zeros(n,K);
        dg=((1:n).^dta.dec)';
        lbnd=0;
        for k=1:dta.K
            ubnd=min(lbnd+dn,dta.n);
            dta.dgsT(lbnd+1:ubnd,k)=FCTR*dg(lbnd+1:ubnd);
            xMax(lbnd+1)=1/dgsT(lbnd+1);
            dta.sT{k}=zeros(dta.n,dta.n);
            for i=lbnd+1:ubnd
                dta.sT{k}(i,i)=FCTR*dg(i);
            end
            lbnd=ubnd;
        end
    else
        dta.sT=cell(dta.K,1);
        lbnd=0;
        for k=1:dta.K
            ubnd=min(lbnd+dn,dta.n);
            ii=lbnd+1:1:ubnd;
            vv=FCTR*ii.^dta.dec;
            xMax(lbnd+1)=1/vv(1);
            dta.sT{k}=sparse(ii,ii,vv,dta.n,dta.n);
            lbnd=ubnd;
        end
    end
    One=0;
    for k=1:K
        One=max(One,norm(dta.sT{k}*xMax)^2);
    end
    dta.One=One;
    dta.NormMax=norm(xMax);
else % "smooth" X 
    dta.box=0;
    dta.sT=cell(K,1);
    dh=2*pi/(n-1);
    TT=zeros(n,n);
    TT(1,1)=1;
    TT(2,1:2)=[-1,1]/dh;
    for i=3:n
        TT(i,[i-2,i-1,i])=[1,-2,1]/dh^2;
    end
    lbnd=0;
    for k=1:K
        dta.sT{k}=sparse([],[],[],n,n);
    end
    for k=1:K
        ubnd=min(lbnd+dn,n);
        if ubnd~=lbnd
            Tmp=zeros(n,n);
            Tmp(lbnd+1:ubnd,:)=TT(lbnd+1:ubnd,:)/sqrt(ubnd-lbnd);
            [ii,jj,vv]=find(Tmp);
            dta.sT{k}=sparse(ii,jj,vv,n,n);
            lbnd=ubnd;
        end
    end
    xx=2*pi*(0:1:n-1)'/(n-1);
    f=1+xx+0.5*xx.^2;
    umx=0;
    for k=1:K
        umx=max(umx,norm(dta.sT{k}*f));
    end
    ff=f/umx;
    dta.NormMax=norm(ff);
    dta.One=0;
    for k=1:K
        dta.One=max(dta.One,norm(dta.sT{k}*ff)^2);
    end
end
end 
%endof create_sT
%

function dta=create_AB(dta,seed)
rng(seed)
if strncmpi(dta.Atype,'g',1)
    dta.Atype='Gaussian';
    dta.A=randn(dta.n);
elseif strncmpi(dta.Atype,'r',1)
    dta.Atype='Rademacher';
    dta.A=sign(randn(dta.n));
else
    dta.Atype='Other';
    [U,~,V]=svd(randn(dta.n));
    d=(1:1:n).^(-dta.edr);
    dta.A=U*diag(d)*V';
end
if strncmpi(dta.Btype,'D',1)
    dta.Btype='Derivative';
    dta.nu=dta.n-1;
    dta.B=zeros(dta.n-1,dta.n);
    for i=1:dta.n-1
        dta.B(i,i)=1;
        dta.B(i,i+1)=-1;
    end
elseif strncmpi(dta.Btype,'P',1)
    dta.Btype='Primitive';
    dta.nu=n;
    dta.B=zeros(n,n);
    for i=1:dta.n
        dta.B(i,1:i)=1;
    end
else
    dta.Btype='Signal';
    dta.nu=dta.n;
    dta.B=eye(dta.n);
end
if ~isfield(dta,'eps')
    dta.eps=0.05;
end
dta.cnr=1/(dta.sigma*icdf('normal',1-dta.eps/2,0,1));
dta.cff=1/dta.cnr^2;
%
dta.Ai=pinv(dta.A);
dta.BAi=dta.B*dta.Ai;
%if dta.box==0,
dta.sTAi=cell(dta.K,1);
for k=1:dta.K
    dta.sTAi{k}=dta.sT{k}*dta.Ai;
end

end


