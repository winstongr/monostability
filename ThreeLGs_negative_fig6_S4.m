clear
set(groot,'DefaultAxesBox','on') %set default, show figure box/frame show
set(groot,'DefaultAxesLinewidth',1) %axis line width
set(groot,'DefaultAxesColor','none') %transparent background
set(groot,'DefaultAxesTicklength',[0.018 0.025])
set(groot,'DefaultFigurePosition',[360,198,560,420]) %figure size, left, bottom, width,height

% to speed up, not using symbolic math
%% initialize and value assignment
Nt=120;                 %Number of Rt values
Rtv=logspace(-1,2,Nt);
Np=42;                  %Number of Cp values
Cpspan=logspace(-1,log10(190),Np);      %Cp range 0.1-190

Nf=35;                  %Number of f values
fv=logspace(0,2,Nf);

D=2;
Kv=[1/3, 0.5, 2];       % DNA binding affinity K*
Nk=length(Kv);
h=2;                    %DNA binding Hill coefficient

Kfn=[1.5,1;1.5,0.25;1.5,0.1;3,0.1];  % for auto-repression, col 1 KN; col 2, fN
Nn=size(Kfn,1);

%assign variables to store results
maxLG=zeros(Np,Nf,Nn,Nk);
maxI=zeros(Np,Nf,Nn,Nk);
LG3_store=zeros(Nt,Np,Nf,Nn,Nk);
LGf_store=zeros(Nt,Np,Nf,Nn,Nk);
LGn_store=zeros(Nt,Np,Nf,Nn,Nk);
LGa_store=zeros(Nt,Np,Nk);
LGb_store=zeros(Nt,Np,Nk);
Rp_store=zeros(Nt,Np,Nk);
Rf_store=zeros(Nt,Np,Nk);

%% calculate Rp, Rf, LG, DNAamount direct not ratio

for ik=1:Nk
    K=Kv(ik);
        for ip=1:Np
            Cp=Cpspan(ip);
            Ct=3;
            for i=1:Nt
                %phos. module
                Rt=Rtv(i);
                Rp=0.5*((Cp+Ct+Rt)-((Cp+Ct+Rt)^2-4*Cp*Rt)^0.5);
                Rp_store(i,ip,ik)=Rp;
                
                LGa=-(Rt*((2*Ct-2*Cp+2*Rt)/(4*((Cp+Ct+Rt)^2-4*Cp*Rt)^(1/2))-1/2))/(Cp/2+Ct/2+Rt/2-((Cp+Ct+Rt)^2-4*Cp*Rt)^(1/2)/2);
                LGa_store(i,ip,ik)=LGa;
                % binding module, solve Rf
                if D==0
                    LGb=1;
                    Rf=Rp;
                    Rf_store(i,ip,ik)=Rf;
                else
                fsol=roots([1  2*D-Rp  K^2  -Rp*K^2]);
                pf=fsol(imag(fsol)==0 & fsol>=0 & fsol<=Rp); % select real, positive Rf solutions
                if length(pf)==1
                    Rf=pf;                                  
                    Rf_store(i,ip,ik)=Rf;
                else
                    Rf=max(pf);
                    Rf_store(i,ip,ik)=-1*length(pf);     %report aberrant fsol (root) with multiple solutions  
                end
                LGb=((K^2+Rf^2)*(K^2+Rf^2+2*D*Rf))/(K^4+2*K^2*Rf^2+4*D*K^2*Rf+Rf^4);
                end
                LGb_store(i,ip,ik)=LGb;
                % autoregulation module
                for j=1:Nf
                    f=fv(j);
                    for in=1:Nn
                        Kn=Kfn(in,1);
                        Rn=Rf/Kn;
                        fn=Kfn(in,2);
                        LGf=(Rf^h*h*(f - 1))/((1 + Rf^h*f)*(1 + Rf^h));
                        LGf_store(i,ip,j,in,ik)=LGf;
                        LGn=(Rn^h*h*(fn - 1))/((1 + Rn^h*fn)*(1 + Rn^h));
                        LGn_store(i,ip,j,in,ik)=LGn;
                        LG3_store(i,ip,j,in,ik)=(LGf+LGn)*LGa*LGb;
                    end
                end
            end
            for j=1:Nf
                for in=1:Nn
                    [Mrt,Irt]=max(LG3_store(:,ip,j,in,ik));
                    maxLG(ip,j,in,ik)=Mrt;
                    maxI(ip,j,in,ik)=Irt;
                end
            end
        end
end


%% Plot LG vs Rt
fLG1=figure('name','LG','DefaultAxesFontSize',14);
fLG1.Position=[360,198,540,250];
%Cpspan(31), 25; fv(13), 5; Kv(2), 0.5; Kfn(1)(3), no CNF or 1.5/0.1
plot(Rtv,LGa_store(:,31,2),'-','linewidth',1,'DisplayName','LGa')
hold on
plot(Rtv,LGb_store(:,31,2),'-','linewidth',1,'DisplayName','LGb')
hold on
plot(Rtv,LGf_store(:,31,13,1,2),'-','linewidth',1,'DisplayName','LGf')
hold on
plot(Rtv,LGn_store(:,31,13,2,2),'-','linewidth',1,'DisplayName','LGn')
hold on
plot([0.1 100],[1 1],'k-','linewidth',1)
hold on
set(gca,'XScale','log')
xlim([0.1 100])
xlabel('Rt')
ylabel('LGs')

fLG2=figure('name','LG','DefaultAxesFontSize',14);
fLG2.Position=[360,198,540,250];
%Cpspan(31), 25; fv(13), 5; Kv(2), 0.5; Kfn(1)(3), no CNF or 1.5/0.1
plot(Rtv,LG3_store(:,31,13,1,2),'k-','linewidth',1,'DisplayName','PF only')
hold on
plot(Rtv,LG3_store(:,31,13,2,2),'b-','linewidth',1,'DisplayName','PF only')
hold on
plot([0.1 100],[1 1],'k-','linewidth',1)
hold on
set(gca,'XScale','log')
xlim([0.1 100])
xlabel('Rt')
ylabel('LG3')
%% plot contour figure

fD=gobjects(Nk,1);
cname=jet(Nn);

for i=1:Nk
    fD(i)=figure('name',append('K',num2str(Kv(i))),'DefaultAxesFontSize',14);
    fD(i).Position=[360,198,540,420];
    for m=1:Nn
        figure(fD(i))
        llabel=['Kn ',num2str(Kfn(m,1)),' fn ',num2str(Kfn(m,2))];
        contour(fv,Cpspan,maxLG(:,:,m,i),[1 1],'-','color',cname(m,:),'linewidth',1,'DisplayName',llabel)
        hold on
    end
    figure(fD(i))
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    title(append('maxLGs K',num2str(Kv(i))))
    ylabel('Cp')
    ylim([0.1 100])
    xlabel('f')
    legend
    colorbar
end

%% plot effect of LGn vs Rp
fLG=figure('name','LG','DefaultAxesFontSize',14);
fLG.Position=[360,198,500,450];
%Cpspan(31), 25; fv(13), 5; Kv(2), 0.5; Kfn(1)(3), no CNF or 1.5/0.1
plot(Rp_store(:,31,2),LG3_store(:,31,13,1,2),'k-','linewidth',1,'DisplayName','PF only')
hold on
plot(Rp_store(:,31,2),LG3_store(:,31,13,2,2),'b-','linewidth',1,'DisplayName','PF only')
hold on
set(gca,'XScale','log')
xlim([0.09 25])
