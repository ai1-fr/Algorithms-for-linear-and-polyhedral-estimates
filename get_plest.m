function res=get_plest(dta,control)
% RES=GET_LPEST(DATA,CONTROL)
% Run NERMLT algorithm to compute polyhedral and linear estimates and dereference the results
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

K=dta.K;
n=dta.n;
nu=dta.nu;
tstart=cputime;
%
res=NERMLTi(dta,control);
res.cpu=cputime-tstart;
%
% dereference estimation matrices
% dereference H matrix of the polyhedral estimate
res.pol.risk=2*sqrt(res.upb);
res.pol.lambda=sqrt(res.upb);
res.pol.mu=res.x/sqrt(res.upb);
res.G=dta.B'*dta.B/res.pol.lambda;
M=zeros(size(res.G));
if dta.box
    dg=(dta.dgsT.^2)*res.pol.mu;
    M=diag(dg);
else
    for k=1:K
        M=M+res.pol.mu(k)*dta.sT{k}'*dta.sT{k};
    end
end
M=0.5*(M+M');
res.G=res.G-M;
res.G=dta.Ai'*res.G*dta.Ai;
res.G=0.5*(res.G+res.G');
[res.U,res.D]=eig(res.G);
res.pol.Theta=1.e-8*eye(dta.n)+res.U*diag(max(diag(res.D),0))*res.U';
res.pol.Theta=0.5*(res.pol.Theta+res.pol.Theta');
[~,P]=eig(res.pol.Theta);
% sTheta=res.U*diag(sqrt(max(diag(P),0)))*res.U';
res.pol.check.Theta=max(0,-min(diag(P)));
Tmp=[res.pol.lambda*eye(dta.nu),0.5*dta.B;0.5*dta.B',dta.A'*res.pol.Theta*dta.A+M];
Tmp=0.5*(Tmp+Tmp');
res.pol.check.LMI=max(0,-min(eig(Tmp)));
res.pol.check.obj=abs(res.pol.lambda+norm(res.pol.mu,dta.q)+dta.cff*trace(res.pol.Theta)-res.pol.risk);
[res.pol.H,res.D]=eig(res.pol.Theta);
res.pol.H=res.pol.H*dta.cnr;
%
sLi=sqrt(1/res.pol.lambda)*eye(nu);
Tmp=zeros(n,n);
for k=1:K
    Tmp=Tmp+res.pol.mu(k)*dta.sT{k}'*dta.sT{k};
end
% dereference H matrix of linear estimate 
Tmp=0.5*(Tmp+Tmp');
[V,D]=eig(Tmp);
sMi=V*diag(1./sqrt(max(diag(D),1.e-16)))*V';
[V,D]=eig(res.pol.Theta);
sTheta=V*diag(sqrt(max(diag(D),0)))*V';
Tmp=sTheta*dta.A*sMi;
[W,D,V]=svd(Tmp);
U=W*V';
S=V*D*V';
Q=sLi*dta.B*sMi*(S+eye(n))^(-1);
res.lin.H=(sqrt(res.pol.lambda)*Q*U'*sTheta)';
F=sqrt(0.5)*(dta.B-res.lin.H'*dta.A)*sMi;
res.lin.mu=res.pol.mu/2;
res.lin.lambda=norm(F)^2;
res.lin.risk=res.lin.lambda+norm(res.lin.mu,dta.q)+dta.sigma*norm(res.lin.H,'fro');
%
if control.print<10000
    fprintf(' %7.1f risks: pol %7.6f lin %7.6f\n',res.cpu,res.pol.risk,res.lin.risk);
    fprintf('This should be near-zero: Theta %.1e LMI %.1e Obj %.1e\n',...
        res.pol.check.Theta,res.pol.check.LMI,res.pol.check.obj);
end

end

