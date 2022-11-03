clear
set(groot,'DefaultAxesBox','on') %set default, show figure box/frame show
set(groot,'DefaultAxesLinewidth',1) %axis line width
set(groot,'DefaultAxesColor','none') %transparent background
set(groot,'DefaultAxesTicklength',[0.018 0.025])
%% initialize conditions


Den=[0.001, 0.75, 1.5, 3, 3.75, 4, 6]; % number of cometing DNA sites, use 0.001 to approximate 0
Nd=length(Den);
Ken=[0.5 2]; % affinity of competing DNA sites
Nk=length(Ken);
h=2;
hl=2;  
f=[5 25];
Nf=length(f);
Rb=[0.75 2.5];
Nb=length(Rb);
NCp=301;        % number of Cp values to scan
Cpv=logspace(-1,3,NCp)';  % Cp values
Ct=3;


Rt_value=zeros(NCp,3,Nd,Nk,Nf,Nb);
RRf_value=zeros(NCp,3,Nd,Nk,Nf,Nb);
Nsol=zeros(NCp,Nd,Nk,Nf,Nb);

%% calulate steady state 
% calculate steady state solutions at different Cp
for n=1:Nb
for m=1:Nf
for j=1:Nk
for i=1:Nd        
    display(['Rb= ',num2str(Rb(n)),' f= ',num2str(f(m)),' Ken= ',num2str(Ken(j)),' and Den= ',num2str(Den(i))]); 
    const=[Den(i) Ken(j) Ct f(m) Rb(n)];
    sol=steady(Cpv,const);
    %i=1;
    Nsol(:,i,j,m,n)=sol(:,1);
    Rt_value(:,:,i,j,m,n)=sol(:,2:4);
    RRf_value(:,:,i,j,m,n)=sol(:,5:7);
    
end
end
end
end
%% calculate occupancy
Krep=[Ken, 1];
Poccup=cell(size(Krep));
for i=1:length(Krep)
    Poccup{i}=RRf_value.^h./(RRf_value.^h+Krep(i)^h);
end

%% figure 7C Occupancy of promoter Krep, Rb=2.5, f=5, K=0.5, different D, Krep (Ken, Kauto)
%(NCp,3,Nd,Nk,Nf,Nb)

clist=jet(Nd);
fig7c=cell(size(Krep));
for m=1:length(Krep)
    occup=Poccup{m};
    fig7c{m}=figure('DefaultAxesFontSize',14); 
    fig7c{m}.Position=[360,198,540,420];
    for i=1:Nd
        occuplot(fig7c{m},Cpv,occup(:,:,i,1,1,2),clist(i,:))
    end
    xlim([0.4 400])
    ylim([-0.1 1.15])
    set(gca,'XScale','log')
    title(['Krep= ',num2str(Krep(m)),' K=0.5 f=5']);

    ylabel('Occupany')
    xlabel('Cp')
end
%% figure 7B  D=4; Rf & occup two Ken two Rb f=5, D=4 match 7A
%(NCp,3,Nd,Nk,Nf,Nb)
% Rf fig 
fig7b=figure('DefaultAxesFontSize',14); 
fig7b.Position=[360,198,220,210];
%fig7b.Position=[360,198,540,450];
clist=["#0072BD","#4DBEEE";"#A2142F","#EDB120"];
di=6; %index for Den
% plot Rf, 
for j=1:Nb
for i=1:Nk
    occuplot(fig7b,Cpv,RRf_value(:,:,di,i,1,j),clist(i,j))
end
end

xlim([0.4 400])
xticks([1 10 100])
ylim([-0.3 5.6])
set(gca,'XScale','log')
title(['fig 7B: D*= ',num2str(Den(di))]);
ylabel('Rf')
xlabel('Cp')

% occupancy fig 
fig7bo=figure('DefaultAxesFontSize',14); 
fig7bo.Position=[360,198,540,450];
occup=Poccup{3};
for j=1:Nb
for i=1:Nk
    occuplot(fig7bo,Cpv,occup(:,:,di,i,1,j),clist(i,j))
end
end

xlim([0.2 400])
ylim([-0.15 1.19])
yticks([0 0.5 1])
set(gca,'XScale','log')

title(['fig 7B: D*= ',num2str(Den(di))]);
ylabel('Occupancy')
xlabel('Cp')
%% fig S5B, Rf, two K D=4, Rb=.75, 2.5
fs5b=figure('DefaultAxesFontSize',14); 
fs5b.Position=[360,198,540,420];
clist=["#0072BD","#4DBEEE";"#A2142F","#EDB120"];
di=6;
for j=1:Nb
for i=1:Nf
    occuplot(fs5b,Cpv,RRf_value(:,:,di,i,1,j),clist(i,j))
end
end

xlim([0.8 400])
ylim([0.05 15])
set(gca,'XScale','log')
set(gca,'YScale','log')
title(['fig S5b: D*= ',num2str(Den(di))]);
ylabel('Rf')
xlabel('Cp')

%% fig S5b side graph, occupancy vs Rf
fs5bo=figure('DefaultAxesFontSize',14); 
fs5bo.Position=[360,198,175,420];
Rfo=logspace(-2,2,101);
occu=cell(length(Krep));

for i=1:length(Krep)
 occu{i}=Rfo.^h./(Krep(i)^h+Rfo.^h); % occupancy of promoter Krep
 plot(occu{i},Rfo,'k-','linewidth',1)
 hold on
end
xlim([-0.16 1.18])
xticks([0 0.5 1])
ylim([0.05 15])
set(gca,'YScale','log')
xlabel('Rf')
ylabel('Occupancy')
title('occu.vs Rf')


%% fig s5C,D,Krep 0.5, 2 , two Rb
fs5=cell(size(Krep));
clist=["#a97c50","#ed1c24";"#3b78bd","#782aff"];
di=6;
for m=1:length(Krep)
    fs5{m}=figure('DefaultAxesFontSize',14); 
    fs5{m}.Position=[360,198,540,420];
    occup=Poccup{m};
    for j=1:Nb
    for i=1:Nk
        occuplot(fs5{m},Cpv,occup(:,:,di,i,1,j),clist(i,j))
    end
    end
    xlim([0.4 400])
    ylim([-0.15 1.19])
    set(gca,'XScale','log')
    yticks([0 0.5 1])
    title(['fig S5CD: Krep= ',num2str(Krep(m))]);
    ylabel('Occupancy')
    xlabel('Cp')
end

%%
function occuplot(x,CpConst,OccuConst,colorN)
Cpv=CpConst;
NCp=length(Cpv);
ssOccup=OccuConst;
bistable=find(ssOccup(:,2)>0,1);
figure(x);
if isempty(bistable)
   plot(Cpv,ssOccup(:,1),'-','color',colorN,'linewidth',1)
   hold on
else
   Ilow=find(ssOccup(:,2)>0,1,'first');
   Ihi=find(ssOccup(:,2)>0,1,'last');
   plot(Cpv(1:Ihi),ssOccup(1:Ihi,1),'-','color',colorN,'linewidth',1)
   hold on
   plot(Cpv(Ilow:Ihi),ssOccup(Ilow:Ihi,2),'-','color',colorN,'linewidth',1)
   hold on
   plot(Cpv(Ilow:Ihi),ssOccup(Ilow:Ihi,3),'-','color',colorN,'linewidth',1)
   hold on
   if Ihi<NCp
      plot(Cpv(Ihi+1:NCp),ssOccup(Ihi+1:NCp,1),'-','color',colorN,'linewidth',1)
      hold on
   end
end

end

function sss=steady(x,const)
Den=const(1);
Ken=const(2);
Ct=const(3);
f=const(4);
Rb=const(5);
n=length(x);
ss=zeros(n,3)-2;
ssp=zeros(n,3);
Nsolution=zeros(n,1);

%h=2;
hl=2;
DNA=Den; %keep the option for h of target promoter
K=Ken;

for i=1:n
    solcount=0;
    Cp=x(i);
    syms RRf Rt
    eqn0=Rb*(1+f*RRf^hl)/(1+(RRf)^hl)-Rt==0;
    if DNA==0
        complx=0;
    else
        complx=sum(RRf^2./(RRf^2+K.^2).*DNA)*2; %RRf amount in DNA complex
    end
    eqn1=0.25*((Cp+Ct+Rt)^2-4*Cp*Rt)==(0.5*(Cp+Ct+Rt)-(RRf+complx))^2;
    eqns=[eqn0 eqn1];
    [Rtv,RRfv]=vpasolve(eqns,[Rt,RRf],[0 inf; 0 Cp],'Random',true);
    
    if (isempty(Rtv))
        Nsolution(i,1)=-1;
    elseif (length(Rtv)==length(RRfv))
        Rtv=sort(Rtv);
        RRfv=sort(RRfv);
        RRp0=0.5*(Cp+Ct+Rtv)-0.5*((Cp+Ct+Rtv).^2-4.*Cp.*Rtv).^0.5;
        for j=1:length(RRfv)
            if DNA==0
                complxv=0;
            else
            complxv=double(subs(complx,RRf,RRfv(j)));
            end
            if (complxv<=RRp0(j))&&(RRfv(j)<=RRp0(j))&&((RRfv(j)+complxv)<=(RRp0(j)+0.1))
                solcount=solcount+1;
                ss(i,solcount)=Rtv(j);
                ssp(i,solcount)=RRfv(j);
            end
        end
        if (solcount==0)
            Nsolution(i,1)=-1;
        else
            Nsolution(i,1)=solcount;
        end
    end
end
sss=[Nsolution ss ssp];   
end

