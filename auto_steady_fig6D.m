clear
set(groot,'DefaultAxesBox','on') %set default, show figure box/frame show
set(groot,'DefaultAxesLinewidth',1) %axis line width
set(groot,'DefaultAxesColor','none') %transparent background
set(groot,'DefaultAxesTicklength',[0.018 0.025])
set(groot,'DefaultFigurePosition',[360,198,560,420])
%% initialize conditions


Den=[0.01, 0.51, 1, 1.51, 1.91, 2.01, 2.11, 2.21, 2.5]; % number of cometing DNA sites, use 0.01 to approximate 0
Nd=length(Den);
Ken=0.5; % affinity of competing DNA sites
h=2;
hl=2;  
f=5;
%Rb=[0.2 1];
Rb=1.52;
NCp=301;        % number of Cp values to scan
Cpv=logspace(0,log10(500),NCp)';  % Cp values
Ct=3;


Rt_value=zeros(NCp,3,Nd);
RRf_value=zeros(NCp,3,Nd);
Nsol=zeros(NCp,Nd);

%% calulate steady state 
% calculate steady state solutions at different Cp
for i=1:Nd        
    const=[Den(i) Ken Ct f Rb];
    sol=steady(Cpv,const);
    %i=1;
    Nsol(:,i)=sol(:,1);
    Rt_value(:,:,i)=sol(:,2:4);
    RRf_value(:,:,i)=sol(:,5:7);
     
end
%% figure 

fssCp=figure('name','SS Rt','DefaultAxesFontSize',14); 
fssCp.Position=[360,198,540,430];
clist=jet(Nd);

for i=1:Nd
    ssRt=Rt_value(:,:,i);
    bistable=find(ssRt(:,2)>0);
    if isempty(bistable)
        plot(Cpv,ssRt(:,1),'-','color',clist(i,:),'linewidth',1)
        hold on
    else
    Ilow=find(ssRt(:,2)>0,1,'first');
    Ihi=find(ssRt(:,2)>0,1,'last');
    plot(Cpv(1:Ihi),ssRt(1:Ihi,1),'-','color',clist(i,:),'linewidth',1)
    hold on
    plot(Cpv(Ilow:Ihi),ssRt(Ilow:Ihi,2),':','color',clist(i,:),'linewidth',1)
    hold on
    plot(Cpv(Ilow:Ihi),ssRt(Ilow:Ihi,3),'-','color',clist(i,:),'linewidth',1)
    hold on
    if Ihi<NCp
        plot(Cpv(Ihi+1:NCp),ssRt(Ihi+1:NCp,1),'-','color',clist(i,:),'linewidth',1)
        hold on
    end
    end
end

xlim([Cpv(1) Cpv(NCp)])
ylim([1 9])
set(gca,'XScale','log')

ylabel('Rt')
xlabel('Cp')


%%



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
