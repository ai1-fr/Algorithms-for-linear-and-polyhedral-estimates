function res=NERMLTi(dta,control)
% Run CTL bundle algorithm
%   By A. Nemirovski 06/27/2023
%

PROX=1;
CHECK=0;
FLAGQ=1;
SMALL=0;
dta.I=control.I;
breaks=[2;1.5;1.1]*control.tolbase;
counter=0;
tstart=cputime;
if strncmpi(control.sol,'m',1)
    cvx_solver mosek
else
    cvx_solver sdpt3
end
cvx_quiet(true)
S=control.S;
K=control.K;
I=control.I;
% n=control.n;
if (control.pi==1)||(K==1)
    pp=2;
else
    lnK=log(K);
    deg=2^floor(log(lnK)/log(2));
    pp=min(2,1+1/deg);
end
itrmax=control.itrmax;
tol=control.tol;
% p=dta.p;
q=dta.q;
R=control.R;
%% minimization over mu: 0\leq mu, \|\mu\|_2\leq R
cfup=0.5;
cflw=0.5;
cflv=0.5;
upb=inf;
lwb=0;
prox=zeros(K,1);
sep.e=zeros(K,1);
sep.b=0;
Bundle=IniBundlePos(K,S,I);
nv=Bundle.nv;
% nc=Bundle.nc;
% dimt=nv-K;
phase=0;
calls=0;
flags=[0;0;0];
res.rep=cell(3,1);
%%
while(1==1)
    phase=phase+1;
    if CHECK
        fprintf('*** phase %d...',phase);
    end
    if phase>1
        if PROX
            prox=xbest;
        end
        as=asbest;
        bs=bsbest;
    end
    x=prox;
    if phase==1
        [as,bs]=GetPiecePos(dta,x);
    end
    calls=calls+1;
    f=dta.cff*sum(max(bs+as*x,0))+norm(x,q);
    %%%%%%%%%%
    if upb>f
        upb=f;
        xbest=x;
        bsbest=bs;
        asbest=as;
    end
    if (control.bmode~='A')||(calls<=2)
        Bundle=AddBundlePos(Bundle,as,bs);
    else
        Bundle=AddBundlePosA(Bundle,as,bs,lmm);
    end
    %
    [~,lin]=omega(prox,pp);
    %%
    cvx_begin
    variable yt(nv,1);
    dual variable lam
    variable obj
    yt >= 0;
    lam: Bundle.Mat*yt <= Bundle.Rhs;
    if FLAGQ
        if q==1
            sum(yt(1:K)) <= R;
        else
            norm(yt(1:K),q) <= R;
        end
    else
        norm(yt(1:K)) <= R;
    end
    if q==1
        sum(yt(1:K))+dta.cff*yt(end) <= obj;
    else
        norm(yt(1:K),q)+dta.cff*yt(end) <= obj;
    end
    obj <= upb;
    minimize obj
    cvx_end
    if CHECK
        fprintf('A: status: %s optval=%7.6f',cvx_status,cvx_optval)
    end
    if ~strncmpi(cvx_status,'S',1)
        fprintf('A: %s',cvx_status);
        if strncmpi(cvx_status,'F',1)
            continue;
        end
    end
    %
    lmm=zeros(Bundle.S,1);
    for s=1:Bundle.S
        lmm(s)=lam((dta.I+1)*s);
    end
    lwb=max(lwb,cvx_optval);
    for fl=1:3
        if (flags(fl)==0)&&(upb<=breaks(fl)*lwb)
            flags(fl)=1;
            res.rep{fl}.phase=phase;
            res.rep{fl}.calls=calls;
            res.rep{fl}.cpu=cputime-tstart;
            res.rep{fl}.x=xbest;
            res.rep{fl}.upb=upb;
            res.rep{fl}.lwb=lwb;
            if control.print<10000
                fprintf('*** %8.1f Tolerance %3.2f phases: %4d calls: %4d lwb: %7.6f upb: %7.6f\n',...
                    res.rep{fl}.cpu,breaks(fl),phase,calls,lwb,upb);
            end
        end
    end
    if (upb<=tol*lwb)|(calls>=control.itrmax)
        res.phase=phase;
        res.calls=calls;
        res.x=xbest;
        res.upb=upb;
        res.lwb=lwb;
        res.cpu=cputime-tstart;
        if control.print<10000,
            fprintf('%8.1f Termination: phases %d calls %d lwb=%7.6f upb=%7.6f\n',...
                res.cpu,phase,calls,lwb,upb);
        end
        return
    end
    %
    lvl=cflv*lwb+(1-cflv)*upb;
    %
    if CHECK
        fprintf('lwb=%7.6f upb=%7.6f lvl=%7.6f',lwb,upb,lvl);
    end
    dup=upb-lvl;
    dlw=lvl-lwb;
    %
    x=prox;
    cvx_begin
    variable yt(nv,1)
    variable dlt
    variable obj
    yt >= 0;
    Bundle.Mat*yt <= Bundle.Rhs;
    if FLAGQ
        if q==1
            sum(yt(1:K)) <= R;
        else
            norm(yt(1:K),q) <= R;
        end
    else
        norm(yt(1:K)) <= R;
    end
    if q==1
        sum(yt(1:K))+dta.cff*yt(end) <= lvl;
    else
        norm(yt(1:K),q)+dta.cff*yt(end) <= lvl;
    end
    
    if pp~=2
        dlt >= norm(yt(1:K),pp);
        minimize 0.5*dlt^2-sum(lin.*yt(1:K));
    else
        minimize norm(prox-yt(1:K));
    end
    cvx_end
    %
    if CHECK
        fprintf('B: status: %s optval=%7.6f',cvx_status,cvx_optval);
    end
    if ~strncmpi(cvx_status,'S',1)
        fprintf('B: %s',cvx_status);
        if ~strncmpi(cvx_status,'F',1)
            continue;
        end
    end
    %
    x=yt(1:K);
    if pp==2
        if cvx_optval <= 1.e-6
            sep.e=zeros(K,1);
            sep.d=0;
        else
            sep.e=yt(1:K)-prox;
            sep.d=norm(sep.e);
            sep.e=sep.e/sep.d;
        end
    else
        if norm(x-prox)<1.e-6
            sep.e=zeros(K,1);
            sep.d=0;
        else
            [~,dw]=omega(x,pp);
            sep.e=dw-lin;
            sep.e=sep.e/norm(sep.e);
            sep.d=sep.e'*(x-prox);
        end
    end
    %
    %
    while(1==1)
        [as,bs]=GetPiecePos(dta,x);
        calls=calls+1;
        f=dta.cff*sum(max(bs+as*x,0))+norm(x,q);
        if CHECK
            G=0.25*dta.B'*dta.B;
            dg=dta.dgsT*x;
            G=G-diag(dg.^2);
            G=dta.Ai'*G*dta.Ai;
            ev=eig(G);
            ff=dta.cff*sum(max(ev,0))+norm(x,q);
            fprintf('f=%7.6f ff=%7.6f',f,ff);
            pause(0.1)
        end
        %
        if upb>f
            upb=f;
            xbest=x;
            asbest=as;
            bsbest=bs;
        end
        for fl=1:3
            if (flags(fl)==0)&(upb<=breaks(fl)*lwb)
                flags(fl)=1;
                res.rep{fl}.phase=phase;
                res.rep{fl}.calls=calls;
                res.rep{fl}.cpu=cputime-tstart;
                res.rep{fl}.x=xbest;
                res.rep{fl}.upb=upb;
                res.rep{fl}.lwb=lwb;
                if control.print<10000
                    fprintf('*** %8.1f Tolerance %3.2f phases: %d calls: %d lwb: %7.6f upb: %7.6f\n',...
                        res.rep{fl}.cpu,breaks(fl),phase,calls,lwb,upb);
                end
            end
        end
        if upb<=tol*lwb
            res.phase=phase;
            res.calls=calls;
            res.x=xbest;
            res.upb=upb;
            res.lwb=lwb;
            res.cpu=cputime-tstart;
            if control.print<10000
                fprintf('%8.1f Termination: phases %d calls %d lwb=%7.6f upb=%7.6f\n',...
                    res.cpu,phase,calls,lwb,upb);
            end
            return;
        end
        if CHECK
            fprintf('calls=%d lwb=%7.6f upb=%7.6f',calls,lwb,upb);
        end
        if control.bmode~='A'
            Bundle=AddBundlePos(Bundle,as,bs);
        else
            Bundle=AddBundlePosA(Bundle,as,bs,lmm);
        end
% 
        if upb-lvl<=cfup*dup
            break;
        end
        if calls>itrmax
            res.phase=phase;
            res.calls=calls;
            res.x=xbest;
            res.upb=upb;
            res.lwb=lwb;
            res.cpu=cputime-tstart;
            if control.print<10000,
                disp(sprintf('%8.1f Termination: phases %d calls %d lwb=%7.6f upb=%7.6f\n',...
                    res.cpu,phase,calls,lwb,upb));
            end
            return;
        end
        %
        cvx_begin
        variable yt(nv,1);
        dual variable lam;
        variable obj;
        yt >= 0;
        lam:Bundle.Mat*yt <= Bundle.Rhs;
        if FLAGQ
            if q==1
                sum(yt(1:K)) <= R;
            else
                norm(yt(1:K),q) <= R;
            end
        else
            norm(yt(1:K)) <= R;
        end
        if q==1
            sum(yt(1:K))+dta.cff*yt(end) <= obj;
        else
            norm(yt(1:K),q)+dta.cff*yt(end) <= obj;
        end
        obj <= lvl;
        sep.e'*(yt(1:K)-prox) >= sep.d-SMALL*abs(sep.d);
        minimize obj;
        cvx_end
        %
        if CHECK
            fprintf('C: status: %s optval=%7.6f',cvx_status,cvx_optval);
        end
        if ~strncmpi(cvx_status,'S',1)
            fprintf('C: %s',cvx_status);
            if ~strncmpi(cvx_status,'I',1)
                break;
            end
            if strncmpi(cvx_status,'F',1)
                continue;
            end
        end
        lmm=zeros(Bundle.S,1);
        for s=1:Bundle.S
            lmm(s)=lam((dta.I+1)*s);
        end
        if cvx_optval>lwb
            lwb=min(cvx_optval,lvl);
            for fl=1:3
                if (flags(fl)==0)&&(upb<=breaks(fl)*lwb)
                    flags(fl)=1;
                    res.rep{fl}.phase=phase;
                    res.rep{fl}.calls=calls;
                    res.rep{fl}.cpu=cputime-tstart;
                    res.rep{fl}.x=xbest;
                    res.rep{fl}.upb=upb;
                    res.rep{fl}.lwb=lwb;
                    if control.print<10000,
                        fprintf('*** %8.1f Tolerance %3.2f phases: %d calls: %d lwb: %7.6f upb: %7.6f\n',...
                            res.rep{fl}.cpu,breaks(fl),phase,calls,lwb,upb);
                    end
                end
            end
            if upb<=tol*lwb
                res.phase=phase;
                res.calls=calls;
                res.x=xbest;
                res.upb=upb;
                res.lwb=lwb;
                res.cpu=cputime-tstart;
                if control.print<10000
                    fprintf('%8.1f Termination: phases %d calls %d lwb=%7.6f upb=%7.6f\n',...
                        res.cpu,phase,calls,lwb,upb);
                end
                return;
            end
            %
            if lvl-lwb<=cflw*dlw
                break;
            end
        end
        %%New point
        cvx_begin
        variable yt(nv,1);
        variable dlt;
        yt >= 0;
        Bundle.Mat*yt <= Bundle.Rhs;
        if FLAGQ
            if q==1
                sum(yt(1:K)) <=R;
            else
                norm(yt(1:K),q) <= R;
            end
        else
            norm(yt(1:K)) <= R;
        end
        if q==1
            sum(yt(1:K))+dta.cff*yt(end) <= lvl;
        else
            norm(yt(1:K),q)+dta.cff*yt(end) <= lvl;
        end
        sep.e'*(yt(1:K)-prox) >= sep.d-SMALL*abs(sep.d);
        if pp~=2
            dlt >= norm(yt(1:K),pp);
            minimize 0.5*dlt^2-sum(lin.*yt(1:K));
        else
            minimize norm(prox-yt(1:K));
        end
        cvx_end
        if CHECK
            fprintf('D: status: %s optval=%7.6f',cvx_status,cvx_optval);
        end
        if ~strncmpi(cvx_status,'S',1)
            fprintf('D: %s',cvx_status);
            if strncmpi(cvx_status,'F',1)
                break;
            end
        end
        x=yt(1:K);
        if pp==2
            if cvx_optval <= 1.e-6
                sep.e=zeros(K,1);
                sep.d=0;
            else
                sep.e=yt(1:K)-prox;
                sep.d=norm(sep.e);
                sep.e=sep.e/sep.d;
            end
        else
            if norm(x-prox) < 1.e-6
                sep.e=zeros(K,1);
                sep.d=0;
            else
                [w,dw]=omega(x,pp);
                sep.e=dw-lin;
                sep.e=sep.e/norm(sep.e);
                sep.d=sep.e'*(x-prox);
            end
        end
        
        
        counter=counter+1;
        if mod(counter,control.print)==0
            tused=cputime-tstart;
            fprintf('%8.1f phases %4d calls %4d lwb=%7.6f upb=%7.6f',...
                tused,phase,calls,lwb,upb);
        end
    end
    if counter&&(mod(counter,control.print)==0)
        tused=cputime-tstart;
        fprintf('%8.1f phases %4d calls %4d lwb=%7.6f upb=%7.6f',...
            tused,phase,calls,lwb,upb);
    end
end
end %endof NERMLT

function Bundle=IniBundlePos(K,S,I);
% s-th element of bundle: \sum_{i=1}^I\max[a_s^ix+b_s^i,0], A^i:1\times K
% Bundle.S: # of elements
% Bundle.K - K
% Bundle.I - I
% Bundle.Mat - constraint matrix of the system
% a_s^ix+b_s^i - t_s^i \leq -b_s^i
% \sum_it_s^i - t \leq 0
% i=1,...,I,s=1,...,S

Bundle.S=S;
Bundle.K=K;
Bundle.I=I;
Bundle.nv=K+I*S+1; %# of variables
Bundle.nc=(I+1)*S; %# of constraints
Bundle.Mat=zeros(Bundle.nc,Bundle.nv);
Bundle.Rhs=zeros(Bundle.nc,1);
row=0;
tbase=K;
for s=1:S
    cbase=row;
    for i=1:I
        row=row+1;
        Bundle.Mat(row,tbase+i)=-1;
        Bundle.Mat(cbase+I+1,tbase+i)=1;
    end
    row=row+1;
    Bundle.Mat(row,end)=-1;
    tbase=tbase+I;
end
Bundle.Sf=0;
Bundle.Rhs=zeros(Bundle.nc,1);
end
% endof IniBundlePos

function [as,bs]=GetPiecePos(dta,mu)
K=dta.K;
I=dta.I;
n=dta.n;
%%%
Mat=dta.BAi'*dta.BAi;
for k=1:K
    Mat=Mat-mu(k)*dta.sTAi{k}'*dta.sTAi{k};
end
Mat=0.5*(Mat+Mat');
[UU,DD]=eig(Mat);
[evs,ind]=sort(-diag(DD));
evs=-evs;
U=zeros(n,n);
D=zeros(n,n);
for i=1:n
    U(:,i)=UU(:,ind(i));
    D(i,i)=D(ind(i),ind(i));
end;
%
As=zeros(n,K);
Bs=zeros(n,1);
BAiU=dta.BAi*U;
AiU=dta.Ai*U;
for i=1:n
    Bs(i)=sum(BAiU(:,i).^2);
end
for k=1:K
    Tmp=dta.sTAi{k}*U;
    for i=1:n
        As(i,k)=-sum(Tmp(:,i).^2);
    end
end
as=zeros(I,K);
bs=zeros(I,1);
for i=1:I
    if i<I
        as(i,:)=As(i,:);
        bs(i)=Bs(i);
    else
        a=zeros(1,K);
        b=0;
        for ii=I:n
            val=As(ii,:)*mu+Bs(ii);
            if val<=0
                continue;
            else
                a=a+As(ii,:);
                b=b+val-As(ii,:)*mu;
            end
        end
        as(i,:)=a;
        bs(i)=b;
    end
end
if 0==1
    zero=inf;
    for i=1:100
        if i==1
            mun=mu;
        else
            10*rand(K,1);
        end
        Matn=dta.BAi'*dta.BAi;
        for k=1:K
            Matn=Matn-mun(k)*dta.sTAi{k}'*dta.sTAi{k};
        end
        Matn=0.5*(Matn+Matn');
        ev=eig(Matn);
        tru=sum(max(ev,0));
        lwb=sum(max(as*mun+bs,0));
        zero=min(zero,tru-lwb);
    end
    if abs(zero)>1.e-6
        zero
    end
end
end %endof GetPiecePos

function Bundle=AddBundlePos(Bnd,as,bs);
% s-th element of bundle: \sum_{i=1}^I\max[a_s^ix+b_s^i], a_s^i: 1\times K
% Bundle.S: # of elements
% Bundle.K - dimension of x
% Bundle.I - I
% Bundle.Mat - constraint matrix of the system
% a_s^ix - t_s^i \leq -b_s^i 
% \sum_it_s^i - t \leq 0
% i=1,...,I,s=1,...,S

S=Bnd.S;
K=Bnd.K;
I=Bnd.I;
Sf=Bnd.Sf;
nv=Bnd.nv;
nc=Bnd.nc;
Bundle=Bnd;
if Sf<Bnd.S
    pos=Sf+1;
    Bundle.Sf=Sf+1;
    cbase=(I+1)*(pos-1);
    Bundle.Sf=pos;
else
    pos=S;
    cbase=(I+1)*(pos-1);
    Bundle.Mat(1:cbase,1:K)=Bundle.Mat(I+2:end,1:K);
    Bundle.Rhs(1:cbase)=Bundle.Rhs(I+2:end);
end
row=cbase;
for i=1:I
    row=row+1;
    Bundle.Mat(row,1:K)=as(i,:);
    Bundle.Rhs(row)=-bs(i);
end
end %endo AddBudlePos

function Bundle=AddBundlePosA(Bnd,as,bs,lmm) % check if this is the good one
% s-th element of bundle: \sum_{i=1}^I\max[a_s^ix+b_s^i], a_s^i: 1\times K
% Bundle.S: # of elements
% Bundle.K - dimension of x
% Bundle.I - I
% Bundle.Mat - constraint matrix of the system
% a_s^ix - t_s^i \leq -b_s^i 
% \sum_it_s^i - t \leq 0
% i=1,...,I,s=1,...,S

S=Bnd.S;
K=Bnd.K;
I=Bnd.I;
Sf=Bnd.Sf;
nv=Bnd.nv;
nc=Bnd.nc;
Bundle=Bnd;
if Sf<Bnd.S
    pos=Sf+1;
    Bundle.Sf=Sf+1;
    cbase=(I+1)*(pos-1);
    Bundle.Sf=pos;
else
    [,ind]=min(lmm);
    pos=ind(1);
    cbase=(I+1)*(pos-1);
end
row=cbase;
for i=1:I
    row=row+1;
    Bundle.Mat(row,1:K+1)=as(i,:);
    Bundle.Rhs(row)=-bs(i);
end
end % endof AddBundlePosA

function [w,dw]=omega(x,pp)
w=0.5*(sum(abs(x).^pp))^(2/pp);
if w==0
    dw=0*x;
else
    dw=(sum(abs(x).^pp)^(2/pp-1))*sign(x).*(abs(x).^(pp-1));
end
end





