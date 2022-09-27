% For contour plots, the script uses the function: othercolor() 
% citation for othercolor:
% Joshua Atkins (2022). othercolor (https://www.mathworks.com/matlabcentral/fileexchange/30564-othercolor), MATLAB Central File Exchange.  

clear
set(groot,'DefaultAxesBox','on') %set default, show figure box/frame show
set(groot,'DefaultAxesLinewidth',1) %axis line width
set(groot,'DefaultAxesColor','none') %transparent background
set(groot,'DefaultAxesTicklength',[0.018 0.025])
set(groot,'DefaultFigurePosition',[360,198,560,420]) %figure size, left, bottom, width,height
%factory default 560,420 for shorter graphs, above 500,420 for thinner graphs
%% initialize, LGa Phos quadratic & LGf autoregulation using symbolic tools
syms Cp Ct Rt
assume(Rt,'positive');
assume(Ct,'positive');
assume(Cp,'positive');
Rp=0.5*((Cp+Ct+Rt)-((Cp+Ct+Rt).^2-4*Cp*Rt).^0.5);
LGa=diff(Rp,Rt)*Rt/(0.5*((Cp+Ct+Rt)-((Cp+Ct+Rt).^2-4*Cp*Rt).^0.5));
dLGa=simplify(diff(LGa,Rt));
%Rtsol=solve(dLGa==0,Rt);

syms Kauto f h Rf
LGf=(Kauto^h*Rf^h*h*(f - 1))/((Kauto^h + Rf^h*f)*(Kauto^h + Rf^h));
sLGf=subs(LGf,[Kauto h],[1 2]);

Nr=51;                  %Number of Rt values for numerical calculation
Rtv=logspace(-1,2,Nr);  % Rt values to scan
cpv=[0.2 1 5];
ncp=length(cpv);
nrf=61;                 % number of Rf values to scan
Rfv=logspace(-1.7,1.3,nrf);

%%  Kauto as unit, cal. Rp/LGa, Cp 0.2,1,5, Ct 1

Rp2_values=zeros(ncp,Nr); %Rp, phosphorylated RR, calculated from Cp, Ct & Rt
LGa2_values=zeros(ncp,Nr);
LGf_value=zeros(nrf,1);

% phosphorylation module
for i=1:ncp
    for j=1:Nr
        Rp2_values(i,j)=double(subs(Rp,[Cp,Ct,Rt],[cpv(i),1,Rtv(j)]));
        LGa2_values(i,j)=double(subs(LGa,[Cp,Ct,Rt],[cpv(i),1,Rtv(j)]));
    end
end
% autoregulation module
for i=1:nrf
    LGf_value(i)=double(subs(sLGf,[f,Rf],[25,Rfv(i)]));
end

%% figure RpRt, LG-Rt, LG-Rp

fRp2=figure('name','RpRt','DefaultAxesFontSize',14);  %Rp-Rt
fLGa2=figure('name','LGaRt','DefaultAxesFontSize',14); %axis Rt
fLGRp=figure('name','LGaRt','DefaultAxesFontSize',14); %axis Rp

for i=1:ncp
    figure(fRp2)
    plot(Rtv,Rp2_values(i,:),'-','linewidth',1,'DisplayName',num2str(cpv(i)))
    hold on
    figure(fLGa2)
    plot(Rtv,LGa2_values(i,:),'-','linewidth',1,'DisplayName',num2str(cpv(i)))
    hold on
    figure(fLGRp)
    plot(Rp2_values(i,:),LGa2_values(i,:),'-','linewidth',1,'DisplayName',num2str(cpv(i)))
    hold on   
end

%% figure touchup
figure(fRp2)
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Rp vs. Rt')
ylabel('Normalized RRp levels (Rp)')
xlabel('Normalized Total RR levels (RT)')
xlim([Rtv(1) Rtv(Nr)])
ylim([0.01 max(cpv)*3])
legend('location','northwest');

figure(fLGa2)
set(gca,'XScale','log')
set(gca,'YScale','linear')
title('LGa-Rt')
ylabel('Logarithmic Gain (LGa)')
xlabel('Normalized Total RR levels (RT)')
xlim([Rtv(1) Rtv(Nr)])
ylim([-0.1 1.2])
legend('location','northwest');

figure(fLGRp)
set(gca,'XScale','log')
set(gca,'YScale','linear')
title('LGa-Rp')
ylabel('Logarithmic Gain (LGa)')
xlabel('Normalized RRp levels (Rp)')
%xlim([0.01 max(cpv)*3])
%ylim([-0.1 1.2])
legend('location','northwest');

%% calculate LGf*LGa, plot LGs
LGLG_values=zeros(ncp,Nr);

for i=1:ncp
    for j=1:Nr
        
        LGLG_values(i,j)=double(subs(sLGf,[f,Rf],[25,Rp2_values(i,j)]))*LGa2_values(i,j);
    end
end

figure(fLGRp)
plot(Rfv',LGf_value,'-','linewidth',1,'DisplayName','LGf')
hold on

for i=1:ncp
    plot(Rp2_values(i,:),LGLG_values(i,:),'-','linewidth',1,'DisplayName',num2str(cpv(i)))
    hold on
end
xlim([0.01 max(cpv)*4])
ylim([-0.1 1.5])

%% Cp=1, Ct varies

Ctv=[0.2 1 5];
Rp_values=zeros(3,Nr);
LGa_values=zeros(3,Nr);

for i=1:length(Ctv)
    for j=1:Nr
        Rp_values(i,j)=double(subs(Rp,[Cp,Ct,Rt],[1,Ctv(i),Rtv(j)]));
        LGa_values(i,j)=double(subs(LGa,[Cp,Ct,Rt],[1,Ctv(i),Rtv(j)]));
    end

end

%% plot figure Ct varies Cp=1
fRp=figure('name','RpRt','DefaultAxesFontSize',14);
fLGa=figure('name','LGaRt','DefaultAxesFontSize',14);

for i=1:length(Ctv)
    figure(fRp)
    plot(Rtv,Rp_values(i,:),'-','linewidth',1,'DisplayName',num2str(Ctv(i)))
    hold on
    figure(fLGa)
    plot(Rtv,LGa_values(i,:),'-','linewidth',1,'DisplayName',num2str(Ctv(i)))
    hold on
end

%% figure Ct touch up
figure(fRp)
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Rp vs. Rt')
ylabel('Normalized RRp levels (Rp/Cp)')
xlabel('Normalized Total RR levels (RT/Cp)')
xlim([Rtv(1) Rtv(Nr)])
ylim([min(Rp_values,[],'all') max(Rp_values,[],'all')*3])
legend('location','northwest');

figure(fLGa)
set(gca,'XScale','log')
set(gca,'YScale','linear')
title('LGa-Rt')
ylabel('Logarithmic Gain (LGa)')
xlabel('Normalized Total RR levels (RT/Cp)')
xlim([Rtv(1) Rtv(Nr)])
ylim([-0.1 1.2])
legend('location','northwest');

%% contour data for Ct Rt
% Rtv 0.1-100, Nr for Rt; Ct 0.1-100, Nr-10
Ctspan=logspace(-1,2,Nr+10);
LGaCt=zeros(Nr,Nr+10);
RpPCt=zeros(Nr,Nr+10);

for i=1:Nr
    for j=1:Nr+10
        RpPCt(i,j)=double(subs(Rp,[Cp,Ct,Rt],[1,Ctspan(j),Rtv(i)]))*100/Rtv(i);
        LGaCt(i,j)=double(subs(LGa,[Cp,Ct,Rt],[1,Ctspan(j),Rtv(i)]));
    end
end

%% contour plot for Ct Rt
fCCt=figure('name','LGa-Ct','DefaultAxesFontSize',14);
fCCt.Position=[360,198,540,420];
mymap=[linspace(0.65,1,64)' linspace(0.65,1,64)' linspace(0.7,1,64)'];
[Mct,cct]=contour(Rtv,Ctspan,LGaCt',[0.5 0.5],'k-','linewidth',1); %can add 'showtext','on' for label 
hold on
pcolor(Rtv,Ctspan,LGaCt');
set(gca,'XScale','log')
set(gca,'YScale','log')
shading interp
colormap(mymap);
colorbar
hold on

title('LGa-Ct')
xlabel('Rt/Cp')
ylabel('Ct/Cp')
%%
fcRp=figure('name','RpPerc','DefaultAxesFontSize',14);
fCRp.Position=[360,198,540,420];

[MRp,cRp]=contour(Rtv,Ctspan,RpPCt',[10 30 50 70 90],'k-','linewidth',1);
hold on
%%
pcolor(Rtv,Ctspan,RpPCt');
set(gca,'XScale','log')
set(gca,'YScale','log')
shading interp
%colormap(flipud(othercolor('RdBu3')))
colormap(flipud(othercolor('RdGy3')))
title('Rp Perc')
ylabel('Ct/Cp')
xlabel('Rt/Cp')
colorbar
%% contour data for Cp-Rt
%Rtv 0.1-100,  Cp 0.1-100
Cpspan=logspace(-1,2,Nr-1);
LGaCp=zeros(Nr,Nr-1);

for i=1:Nr
    for j=1:Nr-1
        LGaCp(i,j)=double(subs(LGa,[Cp,Ct,Rt],[Cpspan(j),1,Rtv(i)]));
    end
end
%% contour plot for Cp Rt
fCCp=figure('name','LGa-Cp','DefaultAxesFontSize',14);
fCCp.Position=[360,198,540,420];
%mymap=[linspace(0.65,1,64)' linspace(0.65,1,64)' linspace(0.7,1,64)'];
[Mcp,ccp]=contour(Rtv,Cpspan,LGaCp',[0.25 0.5],'k-','linewidth',1); %can add 'showtext','on' for label 
hold on

%%
pcolor(Rtv,Cpspan,LGaCp');
set(gca,'XScale','log')
set(gca,'YScale','log')
shading interp
colormap(mymap);
colorbar
hold on

title('LGa-Cp')
ylabel('Cp')
xlabel('Rt')


