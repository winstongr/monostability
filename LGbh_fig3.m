clear
set(groot,'DefaultAxesBox','on') %set default, show figure box/frame show
set(groot,'DefaultAxesLinewidth',1) %axis line width
set(groot,'DefaultAxesColor','none') %transparent background
set(groot,'DefaultAxesTicklength',[0.018 0.025])
set(groot,'DefaultFigurePosition',[360,198,560,420]) %figure size, left, bottom, width,height
%factory default 560,420 for shorter graphs, above 500,420 for thinner graphs
%% initialize for binding module
% DNA amount use direct amount, not ratio of Rt, because there is no Rt
% value if phos. module is not involved
DNAv=[0 0.3 1 3];
Kv=[1/3, 0.5, 1, 3];
nd=length(DNAv);
nk=length(Kv);
n=61;
Rfv=logspace(-1.7,1.3,n)';

Rp_values=zeros(n,nd,nk);
LGb_values=zeros(n,nd,nk);
LGbLGf=zeros(n,nd,nk);

%% calulate Rp, LGb values, LGb*LGf
Rf=Rfv;
f=5;
ha=2; %ha, autoreg. Hill coefficient
LGf=(Rf.^ha*ha*(f - 1))./((1 + Rf.^ha*f).*(1 + Rf.^ha));
h=2; %h, decoy site Hill coefficient
for i=1:nd
    D=DNAv(i);
    for j=1:nk
        K=Kv(j);
        %LGb=((K^2+Rf.^2).*(K^2+Rf.^2+2*D.*Rf))./(K^4+2*K^2.*Rf.^2+4*D*K^2.*Rf+Rf.^4);
        LGb=((K^h+Rf.^h).*(2*D.*Rf.^(h-1)+K^h+Rf.^h))./(K^(2*h)+Rf.^(2*h)+2*K^h.*Rf.^h+2*D*K^h*Rf.^(h-1)*h);  
        LGb_values(:,i,j)=LGb;
        LGbLGf(:,i,j)=LGb.*LGf;
        Rp=Rf+2*D.*Rf.^h./(K^h+Rf.^h);
        Rp_values(:,i,j)=Rp;
    end
end

%% plot S1 Fig:  Rp,LGb, LGb*LGf vs Rf, with varying DNA
fDf=gobjects(nk,1);
fDp=gobjects(nk,1);
cDname=jet(nd);

for i=2:2
    fDf(i)=figure('name',append('LG-Rf K',num2str(Kv(i))),'DefaultAxesFontSize',14);
    fDp(i)=figure('name',append('LG-Rp K',num2str(Kv(i))),'DefaultAxesFontSize',14);
    for j=1:nd
        figure(fDf(i))
        plot(Rfv,LGb_values(:,j,i),'--','linewidth',1,'color',cDname(j,:),'DisplayName',append('D',num2str(DNAv(j))))
        hold on
        plot(Rfv,LGbLGf(:,j,i),'-','linewidth',1,'color',cDname(j,:),'DisplayName',append('D',num2str(DNAv(j))))
        hold on
        figure(fDp(i))
        plot(Rp_values(:,j,i),LGb_values(:,j,i),'--','linewidth',1,'color',cDname(j,:),'DisplayName',append('D',num2str(DNAv(j))))
        hold on
        plot(Rp_values(:,j,i),LGbLGf(:,j,i),'-','linewidth',1,'color',cDname(j,:),'DisplayName',append('D',num2str(DNAv(j))))
        hold on
    end
    figure(fDf(i))
    plot([Rfv(1), Rfv(n)],[1,1],'k:')
    hold on
    set(gca,'XScale','log')
    set(gca,'YScale','linear')
    ylim([-0.07 3])
    xlim([Rfv(1),max(Rfv)])
    title(append('LGs K',num2str(Kv(i))))
    ylabel('LGs')
    xlabel('Rf')
    legend('location','northwest');
    
    figure(fDp(i))
    plot([Rfv(1), Rfv(n)],[1,1],'k:')
    hold on
    set(gca,'XScale','log')
    set(gca,'YScale','linear')
    ylim([-0.07 3])
    xlim([min(Rp_values,[],'all'),max(Rp_values,[],'all')])
    title(append('LGs K',num2str(Kv(i))))
    ylabel('LGs')
    xlabel('Rp')
    legend('location','northwest');
end
%%  plot with vary
fKf=gobjects(nd,1);
fKp=gobjects(nd,1);
cKname=jet(nk);

for i=3:3
    fKf(i)=figure('name',append('LG-Rf D',num2str(DNAv(i))),'DefaultAxesFontSize',14);
    fKp(i)=figure('name',append('LG-Rp D',num2str(DNAv(i))),'DefaultAxesFontSize',14);
    for j=1:nk
        figure(fKf(i))
        plot(Rfv,LGb_values(:,i,j),'--','linewidth',1,'color',cKname(j,:),'DisplayName',append('K',num2str(Kv(j))))
        hold on
        plot(Rfv,LGbLGf(:,i,j),'-','linewidth',1,'color',cKname(j,:),'DisplayName',append('K',num2str(Kv(j))))
        hold on
        figure(fKp(i))
        plot(Rp_values(:,i,j),LGb_values(:,i,j),'--','linewidth',1,'color',cKname(j,:),'DisplayName',append('K',num2str(Kv(j))))
        hold on
        plot(Rp_values(:,i,j),LGbLGf(:,i,j),'-','linewidth',1,'color',cKname(j,:),'DisplayName',append('K',num2str(Kv(j))))
        hold on
    end
    figure(fKf(i))
    plot([Rfv(1), Rfv(n)],[1,1],'k:')
    hold on
    set(gca,'XScale','log')
    set(gca,'YScale','linear')
    ylim([-0.07 3])
    xlim([Rfv(1),max(Rfv)])
    title(append('LGs D',num2str(DNAv(i))))
    ylabel('LGs')
    xlabel('Rf')
    legend('location','northwest');
    
    figure(fKp(i))
    plot([Rfv(1), Rfv(n)],[1,1],'k:')
    hold on
    set(gca,'XScale','log')
    set(gca,'YScale','linear')
    ylim([-0.07 3])
    xlim([min(Rp_values,[],'all'),max(Rp_values,[],'all')])
    title(append('LGs D',num2str(DNAv(i))))
    ylabel('LGs')
    xlabel('Rp')
    legend('location','northwest');
end
%% plot Rp vs Rf, with different DNA or K
fRpD=figure('name','RfRp','DefaultAxesFontSize',14);

for i=1:nd
    figure(fRpD)
    plot(Rp_values(:,i,2),Rfv,'-','linewidth',1,'DisplayName',append('D',num2str(DNAv(i))))
    hold on
end


figure(fRpD)
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Rp vs. Rf')
xlim([0.1,max(Rp_values,[],'all')])
plot([0.1,max(Rp_values,[],'all')],[Kv(2) Kv(2)],'b-','linewidth',1)
ylabel('Rf')
xlabel('Rp')
legend('location','northwest');

%% plot selected LGb and LGb*LGf
fLGsRf=figure('name','LGs-Rf','DefaultAxesFontSize',14);
plot(Rfv,LGf,'k-','linewidth',1,'DisplayName','LGf')
hold on
plot([Rfv(1), Rfv(n)],[1,1],'k--')
hold on

plot(Rfv,LGb_values(:,3,1),'--','linewidth',1,'color','m','DisplayName',append('K',num2str(Kv(1)),'D1'))
hold on
plot(Rfv,LGbLGf(:,3,1),'-','linewidth',1,'color','m','DisplayName',append('K',num2str(Kv(1)),'D1'))
hold on

plot(Rfv,LGb_values(:,4,3),'--','linewidth',1,'color','y','DisplayName',append('K',num2str(Kv(3)),'D3'))
hold on
plot(Rfv,LGbLGf(:,4,3),'-','linewidth',1,'color','y','DisplayName',append('K',num2str(Kv(3)),'D3'))
hold on

set(gca,'XScale','log')
set(gca,'YScale','linear')
ylim([-0.07 3])
xlim([Rfv(1),max(Rfv)])
title('LGs D-K')
ylabel('LGs')
xlabel('Rf')
legend('location','northwest');

%% contour values cal. DNA, f
Rfv2=logspace(-1.7,2,n)';
Nf=41; %f values
fv=logspace(0,2,Nf)';
Nd2=31;
%Kv2=logspace(-1,log10(5),Nk2)';
DNAv2=logspace(-1,1,Nd2);
Kv=[1/3, 0.5, 1, 2];
nk=length(Kv);

maxLGs=zeros(Nd2,Nf,nk);
LG2s=zeros(n,Nd2,Nf,nk);

LGb_store=zeros(n,Nd2,nk);
LGf_store=zeros(n,Nf);

Rf=Rfv2;
ha=2; %ha, autoreg. Hill coefficient
h=2; %h, decoy site Hill coefficient


for i=1:nk
    K=Kv(i);
    for j=1:Nd2
        D=DNAv2(j);
        LGb=((K^h+Rf.^h).*(2*D.*Rf.^(h-1)+K^h+Rf.^h))./(K^(2*h)+Rf.^(2*h)+2*K^h.*Rf.^h+2*D*K^h*Rf.^(h-1)*h);  
        LGb_store(:,j,i)=LGb;
        for m=1:Nf
            f=fv(m);
            LGf=(Rf.^ha*ha*(f - 1))./((1 + Rf.^ha*f).*(1 + Rf.^ha));
            LGf_store(:,m)=LGf;
            LG2=LGb.*LGf;
            LG2s(:,j,m,i)=LG2;
            maxLGs(j,m,i)=max(LG2);
        end
    end
end   

%% contour plot 
fCon=figure('name','LGLG-fDNA','DefaultAxesFontSize',14);
fCon.Position=[360,198,540,420];
%mymap=[linspace(0.65,1,64)' linspace(0.65,1,64)' linspace(0.7,1,64)'];
lspec={'k-','r-','c-','g-','b-'};

pcolor(fv,DNAv2,maxLGs(:,:,4))
set(gca,'XScale','log')
set(gca,'YScale','log')
shading interp
colormap(flipud(othercolor('RdGy3')))
%colormap(mymap)
colorbar
hold on

i=2;
contour(fv,DNAv2,maxLGs(:,:,i),[1 1],lspec{i},'linewidth',1);
hold on
i=4;
contour(fv,DNAv2,maxLGs(:,:,i),[1 1],lspec{i},'linewidth',1);
hold on
    
%set(gca,'XScale','log')
%set(gca,'YScale','log')
title('LGaLGf-DNA&f')
ylabel('DNA')
xlabel('f')

%% to compare h effect
hv=[1.5 2 2.5];
nh=length(hv);

maxLGh=zeros(Nd2,Nf,nh,nk);
LG2h=zeros(n,Nd2,Nf,nh,nk);

LGbh_store=zeros(n,Nd2,nh,nk);
LGf_store=zeros(n,Nf);

Rf=Rfv2;
ha=2; %ha, autoreg. Hill coefficient

for l=1:nk
    K=Kv(l);
for i=1:nh
    h=hv(i);
    for j=1:Nd2
        D=DNAv2(j);
        LGb=((K^h+Rf.^h).*(2*D.*Rf.^(h-1)+K^h+Rf.^h))./(K^(2*h)+Rf.^(2*h)+2*K^h.*Rf.^h+2*D*K^h*Rf.^(h-1)*h);  
        LGbh_store(:,j,i,l)=LGb;
        for m=1:Nf
            f=fv(m);
            LGf=(Rf.^ha*ha*(f - 1))./((1 + Rf.^ha*f).*(1 + Rf.^ha));
            LGf_store(:,m)=LGf;
            LG2=LGb.*LGf;
            LG2h(:,j,m,i,l)=LG2;
            maxLGh(j,m,i,l)=max(LG2);
        end
    end
end 
end

%% contour plot different h
fCon2=figure('name','LGLG-fDNA','DefaultAxesFontSize',14);
fCon2.Position=[360,198,540,420];

lspec={'k-','r-','c-','g-','b-'};
lsp2={'k--','r--','c--','g--','b--'};

pcolor(fv,DNAv2,maxLGh(:,:,2,2))
set(gca,'XScale','log')
set(gca,'YScale','log')
shading interp
colormap(flipud(othercolor('RdBu3')))
colorbar
hold on


for i=1:nh
    contour(fv,DNAv2,maxLGh(:,:,i,2),[1 1],lspec{i},'linewidth',1);
    hold on
    contour(fv,DNAv2,maxLGh(:,:,i,3),[1 1],lsp2{i},'linewidth',1);
    hold on
end
    
%set(gca,'XScale','log')
%set(gca,'YScale','log')
title('LGaLGf-DNA&f')
ylabel('DNA')
xlabel('f')

%% plot LGb-Rf different h
% K=1, D=1 or 2.93 (16 or 23), three h
fLGh=figure('name','LG-Rf-h','DefaultAxesFontSize',14);
for i=1:nh
    plot(Rfv2,LGbh_store(:,16,i,2),lspec{i},'linewidth',1)
    hold on
    plot(Rfv2,LG2h(:,16,15,i,2),lspec{i},'linewidth',1)
    hold on
end
plot(Rfv2,LGf_store(:,15),'k--')
hold on
plot([Rfv2(1), Rfv2(n)],[1,1],'k--')
hold on
set(gca,'XScale','log')
set(gca,'YScale','linear')
ylim([-0.07 2.2])
xlim([Rfv(1),max(Rfv)])
