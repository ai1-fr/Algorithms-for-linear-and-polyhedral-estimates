function res=IPMT(dta,control)
% RES=IPMT(DATA,CONTROL)
% Use IPM to compute polyhedral and linear estimates and dereference the results
%
% Inputs: 
%   DATA: data structure as built in SET_CTL
%   CONTROL: control structure as built in SET_CTL
% Outputs:
%   RES: structure with fields
%   RES.POL: parameters of the computed polyhedral estimate
%   RES.LIN: parameters of the computed linear estimate
%   By A. Nemirovski 06/27/2023
%

if strncmpi(control.sol,'m',1) 
    cvx_solver mosek
else 
    cvx_solver SDPT3
end
if dta.n<=64
    cvx_quiet(true);
else
    cvx_quiet(false);
end

n=dta.n;
nu=dta.nu;
K=dta.K;
% T=zeros(n^2,K);
q=dta.q;
for k=1:K
    if dta.box
        Ts(:,k)=reshape(diag(dta.dgsT(:,k).^2),n^2,1);
    else
        Ts(:,k)=reshape(dta.sT{k}'*dta.sT{k},n^2,1);
    end
end
if 1==1
    tstart=cputime;
    cvx_begin
    variable Theta(n,n) symmetric;
    variable mmu(K,1);
    variable lam;
    variable t;
    mmu >= 0;
    Theta == semidefinite(n);
    [lam*eye(nu),dta.B;dta.B',dta.A'*Theta*dta.A+reshape(Ts*mmu,n,n)] == semidefinite(nu+n);
    t >= lam+norm(mmu,q)+dta.cff*trace(Theta);
    minimize t;
    cvx_end
    res.pol.cpu=cputime-tstart;
    res.pol.status=cvx_status;
    res.pol.risk=cvx_optval;
    res.pol.lambda=lam;
    res.pol.mu=mmu;
    res.pol.Theta=Theta;
    [res.pol.H,D]=eig(Theta);
    res.pol.H=res.pol.H*dta.cnr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart=cputime;
cvx_begin
variable mmu(K,1);
variable lam;
variable H(n,nu);
variable t;
mmu >= 0;
[lam*eye(nu),0.5*(dta.B-H'*dta.A);0.5*(dta.B'-dta.A'*H),reshape(Ts*mmu,n,n)] == semidefinite(nu+n);
t >= lam+norm(mmu,q)+dta.sigma*norm(H,'fro');
minimize t;
cvx_end
res.lin.status=cvx_status;
res.lin.cpu=cputime-tstart;
res.lin.H=H;
res.lin.risk=t;
res.lin.mu=mmu;
if control.print<10000
    fprintf('%4.1f Pol: %s risk %7.6f\n',res.pol.cpu,res.pol.status,res.pol.risk);
    fprintf('%4.1f Lin: %s risk %7.6f\n',res.lin.cpu,res.lin.status,res.lin.risk);
end