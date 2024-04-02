% test CTL algorithm for computing parameters of linear and polyhedral estimates
% Experiments from the paper "First order algorithms for computing linear and polyhedral estimates"
% by Y. Bekri, A. Juditsky, and A. Nemirovski
% Version 3/31/2024

clear dta cntr res_ctl res_ipm
ns=[64 128 256 512 1024 1024];
Ks=[8 16 32 64 128 1024];
ne=length(ns);

calls=zeros(1,ne);
phases=zeros(1,ne);
cpu_ctl=zeros(1,ne);
erisk=zeros(1,ne);

nipm=find(ns<=64);
cpu_ipm=zeros(1,nipm);
eripm=zeros(1,nipm);
for ie=1:ne
    [dta,cntr] = set_ctl(ns(ie),Ks(ie)); % set parameters of the experiment and define the data of the estimation model    
     if ie==1
         fprintf(' Running CTL: n=%4d, K=%4d', ns(ie), Ks(ie))
         pause(0.1)
     else
         fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b Running CTL: n=%4d, K=%4d', ns(ie), Ks(ie))
         pause(0.1)
     end
    res_ctl=get_plest(dta,cntr); % run NERMLT algorithm and dereference results
    cpu_ctl(ie)=res_ctl.cpu;
    calls(ie)=res_ctl.calls;
    phases(ie)=res_ctl.phase;
    erisk(ie)=res_ctl.pol.risk;
    if ns(ie)<=64
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b Running CVX: n=%4d, K=%4d', ns(ie), Ks(ie))
        pause(0.1)
        res_ipm=IPMT(dta,cntr); % run IPM solver by CVX
        cpu_ipm(ie)=res_ipm.pol.cpu;
        eripm(ie)=res_ipm.pol.risk;
    end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\n')

% results for the simulation with 6 experiments 
fprintf(' Results for the CTL algorithm\n')
fprintf('     n    %4d    %4d     %4d     %4d     %4d     %4d\n', ns)
fprintf('     K    %4d    %4d     %4d     %4d     %4d     %4d\n', Ks)
fprintf(' calls    %4d    %4d     %4d     %4d     %4d     %4d\n', calls)
fprintf('phases    %4d    %4d     %4d     %4d     %4d     %4d\n', phases)
fprintf(' CPU s   %5.1f   %5.1f    %5.1f    %5.1f    %5.1f    %5.1f\n', cpu_ctl)
fprintf(' Risks  %5.4f  %5.4f   %5.4f   %5.4f   %5.4f   %5.4f\n', erisk)

fprintf(' Results for the IPM sover by CVX; n=64\n')
fprintf(' CPU s   %5.1f\n', cpu_ipm) 
fprintf(' Risk   %5.4f\n',eripm)






