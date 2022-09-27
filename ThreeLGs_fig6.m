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


DNAamount=[0 0.25, 0.5,1,1.5, 2];   % DNA amount D*
Nd=length(DNAamount); 

Kv=[1/3, 0.5, 2];       % DNA binding affinity K*
Nk=length(Kv);
h=2;                    %DNA binding Hill coefficient

%assign variables to store results
maxLG=zeros(Np,Nf,Nd,Nk);
maxRpPerc=zeros(Np,Nf,Nd,Nk);
maxI=zeros(Np,Nf,Nd,Nk);
LG3_store=zeros(Nt,Np,Nf,Nd,Nk);
LGf_store=zeros(Nt,Np,Nf,Nd,Nk);
LGfb_store=zeros(Nt,Np,Nf,Nd,Nk);
LGa_store=zeros(Nt,Np,Nd,Nk);
LGb_store=zeros(Nt,Np,Nd,Nk);
Rp_store=zeros(Nt,Np,Nd,Nk);
RpPerc_store=zeros(Nt,Np,Nd,Nk);
Rf_store=zeros(Nt,Np,Nd,Nk);

%% calculate Rp, Rf, LG, DNAamount direct not ratio

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
                RpPerc_store(i,ip,id,ik)=Rp*100/Rt;
                LGa=-(Rt*((2*Ct-2*Cp+2*Rt)/(4*((Cp+Ct+Rt)^2-4*Cp*Rt)^(1/2))-1/2))/(Cp/2+Ct/2+Rt/2-((Cp+Ct+Rt)^2-4*Cp*Rt)^(1/2)/2);
                LGa_store(i,ip,id,ik)=LGa;
                % binding module, solve Rf
                if D==0
                    LGb=1;
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
                LGb=((K^2+Rf^2)*(K^2+Rf^2+2*D*Rf))/(K^4+2*K^2*Rf^2+4*D*K^2*Rf+Rf^4);
                end
                LGb_store(i,ip,id,ik)=LGb;
                % autoregulation module
                for j=1:Nf
                    f=fv(j);
                    LGf=(Rf^h*h*(f - 1))/((1 + Rf^h*f)*(1 + Rf^h));
                    LGf_store(i,ip,j,id,ik)=LGf;
                    LGfb_store(i,ip,j,id,ik)=LGf*LGb;
                    LG3_store(i,ip,j,id,ik)=LGf*LGa*LGb;
                end
            end
            for j=1:Nf
                [Mrt,Irt]=max(LG3_store(:,ip,j,id,ik));
                maxLG(ip,j,id,ik)=Mrt;
                maxI(ip,j,id,ik)=Irt;
                maxRpPerc(ip,j,id,ik)=RpPerc_store(Irt,ip,id,ik);
            end
        end
    end
end

%% plot LGs figure
% maxLG vs Cp-f. 1, DNA varies; 2, K varies; 3, selected few DNA/K
% maxRpPerc vs Cp-f, same as above series
fLGa=gobjects(Nk,1);
fLG2=gobjects(Nk,1);
cname=jet(Nd);

for i=2:2
    fLGa(i)=figure('name',append('K',num2str(Kv(i))),'DefaultAxesFontSize',14);
    fLGa(i).Position=[360,198,540,250];
    Cpip=[14,22,31];
    for j=1:length(Cpip)
        figure(fLGa(i))
        plot(Rp_store(:,Cpip(j),1,Nk),LGa_store(:,Cpip(j),1,Nk),'k:','linewidth',1)
        hold on
        plot(Rp_store(:,Cpip(j),Nd,i),LG3_store(:,Cpip(j),13,Nd,i),'-','color',cname(Nd,:),'linewidth',1)
        hold on
    end
    set(gca,'XScale','log')
    set(gca,'YScale','linear')
    xlim([0.04 40])
    ylim([-0.05 1.7])
    title(append('LGa K',num2str(Kv(i))))
    ylabel('LGa')
    xlabel('Rp')
    
    fLG2(i)=figure('name',append('K',num2str(Kv(i))),'DefaultAxesFontSize',14);
    fLG2(i).Position=[360,198,540,250];
    for m=1:Nd
        figure(fLG2(i))
        plot(Rp_store(:,Np,m,i),LGfb_store(:,Np,13,m,i),'--','color',cname(m,:),'linewidth',1)
        hold on
        plot(Rp_store(:,2,m,i),LGfb_store(:,2,13,m,i),'--','color',cname(m,:),'linewidth',1)
        hold on
        %plot(Rp_store(:,Np,m,i),LG3_store(:,Np,20,m,i),'-','color',cname(m,:),'linewidth',1)
        %hold on
    end
    figure(fLG2(i))
    set(gca,'XScale','log')
    set(gca,'YScale','linear')
    title(append('LGs K',num2str(Kv(i))))
    xlim([0.04 40])
    ylim([-0.05 1.7])
    ylabel('LGf*LGb')
    xlabel('Rp')
end


%% plot contour figure

fD=gobjects(Nk,1);
cname=jet(Nd);

for i=2:2
    fD(i)=figure('name',append('K',num2str(Kv(i))),'DefaultAxesFontSize',14);
    fD(i).Position=[360,198,540,420];
    for m=1:Nd
        figure(fD(i))
        contour(fv,Cpspan,maxLG(:,:,m,i),[1 1],'-','color',cname(m,:),'linewidth',1)
        hold on
    end
    figure(fD(i))
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    title(append('maxLGs K',num2str(Kv(i))))
    ylabel('Cp')
    ylim([0.1 100])
    xlabel('f')
    colorbar
end

