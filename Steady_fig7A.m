clear
set(groot,'DefaultAxesBox','on') %set default, show figure box/frame show
set(groot,'DefaultAxesLinewidth',1) %axis line width
set(groot,'DefaultAxesColor','none') %transparent background
set(groot,'DefaultAxesTicklength',[0.018 0.025])
set(groot,'DefaultFigurePosition',[360,198,560,420]) %figure size, left, bottom, width,height

% to speed up, not using symbolic math
%% initialize and value assignment
Nt=120;                 %Number of Rt values
Rtv=logspace(-1,2,Nt)';
Np=38;                  %Number of Cp values
Cpspan=logspace(-1,log10(500),Np);      %Cp range 0.1-500

DNAamount=[0, 4];   % DNA amount D*
Nd=length(DNAamount); 

Kv=[0.5, 2];       % DNA binding affinity K*
Nk=length(Kv);
Rbv=[0.75 2.5];
Nb=length(Rbv);
fv=[5, 25];
Nf=length(fv); 
h=2;                    %DNA binding Hill coefficient

%assign variables to store results
Rp_store=zeros(Nt,Np,Nd,Nk);
Rf_store=zeros(Nt,Np,Nd,Nk);
Frt=zeros(Nt,Np,Nd,Nk,Nf,Nb);

%% calculate Rp, Rf, Occupancy based on Phos./DNA binding modules with known Rt 

for ik=1:Nk
    K=Kv(ik);
    for id=1:Nd
        D=DNAamount(id);
        for ip=1:Np
            Cp=Cpspan(ip);
            Ct=3;
            for i=1:Nt
                %phos. module
                Rt=Rtv(i);
                Rp=0.5*((Cp+Ct+Rt)-((Cp+Ct+Rt)^2-4*Cp*Rt)^0.5);
                Rp_store(i,ip,id,ik)=Rp;
                % binding module, solve Rf
                if D==0
                    Rf=Rp;
                    Rf_store(i,ip,id,ik)=Rf;
                else
                fsol=roots([1  2*D-Rp  K^2  -Rp*K^2]);
                pf=fsol(imag(fsol)==0 & fsol>=0 & fsol<=Rp); % select real, positive Rf solutions
                    if length(pf)==1
                        Rf=pf;                                  
                        Rf_store(i,ip,id,ik)=Rf;
                    else
                        Rf=max(pf);
                        Rf_store(i,ip,id,ik)=-1*length(pf);     %report aberrant fsol (root) with multiple solutions  
                    end
                end
                for ib=1:Nb
                    for ii=1:Nf
                        Frt(i,ip,id,ik,ii,ib)=Rbv(ib)*(1+fv(ii)*Rf^h)/(1+Rf^h);
                    end
                end
            end
            
        end
    end
end
occup=Rf_store.^h./(Rf_store.^h+1);

%% plot figure selected f, K, R
%(Nd,Nk,Nf,Nb);
j=31;
fig7a=figure('name','fig7A','DefaultAxesFontSize',14);
fig7a.Position=[360,198,480,440];
plot(Rtv,Frt(:,j,2,1,1,2),'r-','linewidth',1) %K 0.5 Rb 2.5 Cp 100
hold on
plot(Rtv,Frt(:,1,2,1,1,2),'b-','linewidth',1)  %K 0.5 Rb 2.5 Cp 0.1
hold on
plot(Rtv,Frt(:,j,2,1,1,1),'r-','linewidth',1)  %K 0.5 Rb 0.75 Cp 100
hold on
plot(Rtv,Frt(:,1,2,1,1,1),'b-','linewidth',1)   %K 0.5 Rb 0.75 Cp 0.1
hold on
plot(Rtv,Frt(:,j,2,2,1,2),'-','linewidth',1)  %K 2 Rb 2.5 Cp 100
hold on
plot(Rtv,Frt(:,1,2,2,1,2),'-','linewidth',1)  %K 2 Rb 2.5 Cp 0.1
hold on
plot(Rtv,Rtv,'k-','linewidth',1)
hold on
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([1 14])
ylim([0.65 14])
title(append('f = ',num2str(fv(1))))
ylabel('F(B(A(Rt)))')
xlabel('Rt')

    



