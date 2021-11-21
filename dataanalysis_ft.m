% 
%  saveas from 150309 spectrometer
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all
close all

fs=16;

namefile='Results Table.csv';
sheet=1;
nameRange='A2:A34'
% nameRange='B2:B4'

[num,namelist,raw]=xlsread(namefile,sheet,nameRange);

filenames=namelist;
subsetAll=[];

for now=1:length(filenames)

    filename = [filenames{now} '.Sample.Raw.csv'];
    sheet = 1;
    wlRange='A2:A132';
    xlRange = 'B2:B132';

    subsetA = xlsread(filename, sheet, xlRange);

    if now==1
        wl = xlsread(filename, sheet, wlRange);
    end

    subsetAll=[subsetAll subsetA];
    
end

dx=[ones(1,33)*1]; % [mm]
conccu=[1*ones(1,3) 0*ones(1,3) 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75 0*ones(1,9) 0.75 0.75 0.75 0.5 0.5 0.5 0.25 0.25 0.25]*1;
concni=[0*ones(1,3) 1*ones(1,3) 0*ones(1,9) 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75]*2.2;
sampleid=[1*ones(1,3) 2*ones(1,3) 3*ones(1,3) 4*ones(1,3) 5*ones(1,3) 6*ones(1,3) 7*ones(1,3) 8*ones(1,3) 9*ones(1,3) 10*ones(1,3) 11*ones(1,3)];

A=subsetAll;

uaA=subsetAll*log(10)./repmat(dx,length(wl),1); % [mm-1]

uaAav=zeros(size(uaA,1),1);
conccuav=[];
concniav=[];
for i=1:max(sampleid)
    uaAav(:,i)=mean(uaA(:,sampleid==i),2);
    conccuav(i)=mean(conccu(sampleid==i));
    concniav(i)=mean(concni(sampleid==i));
end

for i=1:max(sampleid)
    uaAmid(:,i)=mean(uaA(:,sampleid==i),2);
    uaAsid(:,i)=std(uaA(:,sampleid==i),[],2);
end
figure;plot(wl,uaAmid)
xlim([750 1100])

figure;plot(wl,uaA)
figure;plot(wl,uaAav)

save('data_sulphates','wl','uaA','uaAav','sampleid','uaAmid','uaAsid','conccuav','concniav')

figure;
subplot(2,1,1)
plot(wl,uaAav(:,conccuav~=0&concniav==0))
title('The copper solutions')
legend('M','0.25 M','0.5 M','0.75 M')
axis tight
ylabel('\mu_a [mm^{-1}]')
subplot(2,1,2)
plot(wl,uaAav(:,conccuav~=0&concniav==0)./repmat(conccuav(conccuav~=0&concniav==0),size(uaA,1),1))
title('The copper solutions normalised')
axis tight
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1} M^{-1}]')

figure;
aux=uaAav(:,conccuav~=0&concniav==0);
plot(wl,[aux(:,2:4) aux(:,1)])
axis tight
ylabel('\mu_a [mm^{-1}]')
l=legend('0.25 M','0.5 M','0.75 M','1 M');
set(l,'linewidth',1','fontsize',fs-5)
axis tight
h=gca;
set(h,'box','off','linewidth',2,'fontsize',fs)
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1}]')


figure;
aux=uaAav(:,conccuav~=0&concniav==0)./repmat(conccuav(conccuav~=0&concniav==0),size(uaA,1),1);
plot(wl,[aux(:,2:4) aux(:,1)])
axis tight
ylabel('\mu_a [mm^{-1}]')
l=legend('0.25 M','0.5 M','0.75 M','1 M');
set(l,'linewidth',1','fontsize',fs-5)
axis tight
h=gca;
set(h,'box','off','linewidth',2,'fontsize',fs)
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1} M^{-1}]')




figure;
subplot(2,1,1)
plot(wl,uaAav(:,conccuav==0&concniav~=0))
title('The nickel solutions')
legend('M','0.25 M','0.5 M','0.75 M')
axis tight
ylabel('\mu_a [mm^{-1}]')
subplot(2,1,2)
plot(wl,uaAav(:,conccuav==0&concniav~=0)./repmat(concniav(conccuav==0&concniav~=0),size(uaA,1),1))
title('The nickel solutions normalised')
xlabel('Wavelength [nm]')
axis tight
ylabel('\mu_a [a.u.]')


figure;
aux=uaAav(:,conccuav==0&concniav~=0);
plot(wl,[aux(:,2:4) aux(:,1)])
axis tight
ylabel('\mu_a [mm^{-1}]')
l=legend('0.55 M','1.1 M','1.65 M','2.2 M');
set(l,'linewidth',1','fontsize',fs-5)
axis tight
h=gca;
set(h,'box','off','linewidth',2,'fontsize',fs)
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1}]')


figure;
aux=uaAav(:,conccuav==0&concniav~=0)./repmat(concniav(conccuav==0&concniav~=0),size(uaA,1),1);
plot(wl,[aux(:,2:4) aux(:,1)])
axis tight
ylabel('\mu_a [mm^{-1}]')
l=legend('0.55 M','1.1 M','1.65 M','2.2 M');
set(l,'linewidth',1','fontsize',fs-5)
axis tight
h=gca;
set(h,'box','off','linewidth',2,'fontsize',fs)
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1} M^{-1}]')



figure;
plot(wl,uaAav(:,conccuav~=0&concniav~=0))
hold on
plot(wl,uaAav(:,1),'color',[0.4 0.4 1],'linewidth',2)
plot(wl,uaAav(:,2),'color',[0.4 1 0.4],'linewidth',2)
title('The mixtures')
legend('25 %','50 %','75%','CuSo4','NiSo4')
axis tight
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1}]')

aux=uaAav(:,conccuav~=0&concniav~=0);

figure;
plot(wl,aux(:,1),'color',[0 0 1],'linewidth',1.5)
hold on
plot(wl,aux(:,2),'color',[0 0.7 0.7],'linewidth',1.5)
plot(wl,aux(:,3),'color',[0 1 0],'linewidth',1.5)
plot(wl,0.75*uaAav(:,1)+0.25*uaAav(:,2),'--','color',0.5*[0 0 1],'linewidth',2)
plot(wl,0.50*uaAav(:,1)+0.50*uaAav(:,2),'--','color',0.5*[0 0.7 0.7],'linewidth',2)
plot(wl,0.25*uaAav(:,1)+0.75*uaAav(:,2),'--','color',0.5*[0 1 0],'linewidth',2)
title('The mixtures')
legend('25 % meas','50 % meas','75 % meas','25 % prev','50 % prev','75 % prev')
axis tight
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1}]')

fs=16;

figure;
plot(wl,uaAav(:,1),'color',[0.6 0.6 1],'linewidth',2)
hold on
plot(wl,uaAav(:,2),'color',[0.4 1 0.4],'linewidth',2)
plot(wl,aux(:,1),'color',[0.2 0.2 1],'linewidth',1.5)
plot(wl,aux(:,2),'color',[0 0.7 0.7],'linewidth',1.5)
plot(wl,aux(:,3),'color',[0 1 0]*0.8,'linewidth',1.5)
plot(wl,0.75*uaAav(:,1)+0.25*uaAav(:,2),'--','color',0.5*[0 0 1],'linewidth',2.5)
plot(wl,0.50*uaAav(:,1)+0.50*uaAav(:,2),'--','color',0.5*[0 0.7 0.7],'linewidth',2.5)
plot(wl,0.25*uaAav(:,1)+0.75*uaAav(:,2),'--','color',0.5*[0 1 0],'linewidth',2.5)
l=legend('0:1, measured','1:0, measured','1:3, measured','2:2, measured','3:1, measured','1:3, predicted','2:2, predicted','1:3, predicted')
set(l,'linewidth',1','fontsize',fs-5)
axis tight
h=gca;
set(h,'box','off','linewidth',2,'fontsize',fs)
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1}]')

%%
% to aid 160511 analysis at wavelengths below 700 nm

Q0=uaAav(:,1)*0.5;
Q100=uaAav(:,2);
figure;
plot(wl,0.75*Q0+0.25*Q100,'--','color',0.5*[0 0 1],'linewidth',2)
hold on
plot(wl,0.50*Q0+0.50*Q100,'--','color',0.5*[0 0.7 0.7],'linewidth',2)
plot(wl,0.25*Q0+0.75*Q100,'--','color',0.5*[0 1 0],'linewidth',2)
plot(wl,Q0,'color',[0.4 0.4 1],'linewidth',2)
plot(wl,Q100,'color',[0.4 1 0.4],'linewidth',2)
title('The mixtures')
legend('hypothetical Q25','hypothetical Q50','hypothetical Q75','Q0=0.5M CuSO_4','Q100=2.2M NiSO_4')
axis tight
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1}]')


% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );
bolorigin=1;
if bolorigin
    opts.Lower = [-Inf 0];
    opts.Upper = [Inf 0];
else
    opts.Lower = [-Inf -Inf];
    opts.Upper = [Inf Inf];
end


for i=1:length(wl)

    aux=uaA(:,conccu==0&concni~=0);
    aux2=concni(conccu==0&concni~=0);
%     aux=[aux(:,2:4) aux(:,1)];
%     aux2=[aux2(2:4) aux2(1)];

    xData=aux2';

    YData=aux;
    yData=YData(i,:)';
    [fitresult, gof] = fit( xData, yData, ft, opts );
    m=coeffvalues(fitresult);
    
    if wl(i)==750
         figure(20);subplot(1,2,1);plot(xData,yData,'o','markersize',7,'markeredgecolor',[0.2 0.5 0.2]);
         hold on;
         plot([0;xData],(m(1)*[0;xData]+m(2)),'--','linewidth',2,'color','r');xlabel('Concentration , M','fontsize',13);ylabel('\mu_a [mm^{-1}]','fontsize',13);title(['R^2:' num2str(gof.rsquare) ', at' ' ' num2str(wl(i)) ' ' 'nm'],'fontsize',11)
         H=gca;set(H,'box','off','linewidth',2,'fontsize',11)
         legend('Data','Linear fit','location','northwest')
         
         figure(200)
         plot(xData,yData,'d','markersize',7,'markerfacecolor',[0.2 0.5 0.2],'markeredgecolor',[0.2 0.5 0.2]);
         hold on;
         plot([0;xData],(m(1)*[0;xData]+m(2)),'--','linewidth',2,'color',[0.3 0.3 0.3]);xlabel('NiSO_4.6H_2O concentration [M]');ylabel(['\mu_{a}(' ' ' num2str(wl(i)) ' ' 'nm)' ' ' '[mm^{-1}]']);
         %title(['R^2:' num2str(gof.rsquare) ', at' ' ' num2str(wl(i)) ' ' 'nm'])
         string_ni=['R^2:' num2str(gof.rsquare)];
         text(1.5,0.4,string_ni,'fontsize',fs-4)
         H=gca;set(H,'box','off','linewidth',2,'fontsize',16)
         l=legend('Data','Linear fit','location','northwest')
         set(l,'fontsize',fs-4)
         axis tight
         
    end
    
    mvec(i,1:2)=m;
    gofvec(i)=gof.rsquare;
    
end

ni_mvec=mvec;
ni_gofvec=gofvec;
    
for i=1:length(wl)

    aux=uaA(:,conccu~=0&concni==0);
    aux2=conccu(conccu~=0&concni==0);
%     aux=[aux(:,2:4) aux(:,1)];
%     aux2=[aux2(2:4) aux2(1)];

    xData=aux2';

    YData=aux;
    yData=YData(i,:)';
    [fitresult, gof] = fit( xData, yData, ft, opts );
    m=coeffvalues(fitresult);
    
    if wl(i)==810
         figure(10);subplot(1,2,1);plot(xData,yData,'o','markersize',7,'markeredgecolor',[0.2 0.2 1]);
         hold on;
         plot([0;xData],(m(1)*[0;xData]+m(2)),'--','linewidth',2,'color','r');xlabel('Concentration , M','fontsize',13);ylabel('\mu_a [mm^{-1}]','fontsize',13);title(['R^2:' num2str(gof.rsquare) ', at' ' ' num2str(wl(i)) ' ' 'nm'],'fontsize',11)
         H=gca;set(H,'box','off','linewidth',2,'fontsize',11)
         legend('Data','Linear fit','location','northwest')
         
         figure(100);
         plot(xData,yData,'d','markersize',7,'markerfacecolor',[0.2 0.2 1],'markeredgecolor',[0.2 0.2 1]);
         hold on;
         plot([0;xData],(m(1)*[0;xData]+m(2)),'--','linewidth',2,'color',[0.3 0.3 0.3]);xlabel('CuSO_4.5H_2O concentration [M]');ylabel(['\mu_{a}(' ' ' num2str(wl(i)) ' ' 'nm)' ' ' '[mm^{-1}]']);
         %title(['R^2:' num2str(gof.rsquare) ', at' ' ' num2str(wl(i)) ' ' 'nm'])
         string_cu=['R^2:' num2str(gof.rsquare)];
         text(0.6,1.2,string_cu,'fontsize',fs-4)
         H=gca;set(H,'box','off','linewidth',2,'fontsize',16)
         l=legend('Data','Linear fit','location','northwest')
         set(l,'fontsize',fs-4)
         axis tight
         
    end
    
    mvec(i,1:2)=m;
    gofvec(i)=gof.rsquare;
    
end

cu_mvec=mvec;
cu_gofvec=gofvec;

figure(10);subplot(1,2,2);
plot(wl,cu_mvec(:,1),'color',[0.2 0.2 1]);
axis tight
h=gca;
set(h,'box','off','linewidth',2)
title('Molar absorption CuSO_4')

figure(20);subplot(1,2,2);
plot(wl,ni_mvec(:,1),'color',[0.2 0.5 0.2]);
axis tight
h=gca;
set(h,'box','off','linewidth',2)
title('Molar absorption NiSO_4')

blood=load('Hbspec.mat');warning('check units')
figure(30);
semilogy(blood.wave,blood.muhhb,'color',[0.5 0 0.5]);
hold on
semilogy(blood.wave,blood.muhbo2,'color',[1 0 0]);
axis tight
h=gca;
set(h,'box','off','linewidth',2)
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1} M^{-1}')
legend('HHb','HbO_2')
title('Blood')
xlim([500 1000])


figure(30);
semilogy(blood.wave,blood.muhhb,'color',[0.5 0 0.5]);
hold on
semilogy(blood.wave,blood.muhbo2,'color',[1 0 0]);
axis tight
h=gca;
set(h,'box','off','linewidth',2)
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1} M^{-1}')
legend('HHb','HbO_2')
title('Blood')
xlim([500 1000])

fs=16;

figure;
[hAx,hLine1,hLine2] = plotyy(wl,cu_mvec(:,1),wl,ni_mvec(:,1))
xlabel('Wavelength [nm]','fontsize',fs+4)
ylabel(hAx(1),'\alpha [mm^{-1} M^{-1}] (CuSO_4.5H_2O)','fontsize',fs+4)
ylabel(hAx(2),'\alpha [mm^{-1} M^{-1}] (NiSO_4.6H_2O)','fontsize',fs+4)
set(hAx,'xlim',[500 1150],'linewidth',2,'box','off','fontsize',fs+4)
set(hAx(1),'ylim',[0 2.7],'ycolor',[0.4 0.4 1],'ytick',[0:0.5:2.5 2.7])
set(hAx(2),'ylim',[0 0.45],'ycolor',[0.2 0.5 0.2]*1.5,'ytick',[0:0.1:0.4 0.45])
set(hLine1,'color',[0.4 0.4 1],'linewidth',2)
set(hLine2,'color',[0.2 0.5 0.2]*1.5,'linewidth',2)
set(gca,'xtick',[500:100:1100])
grid on

% hold on
% plot(wl,ni_mvec(:,1))
% xlim([500 1150])

figure
[hAx3,hLine3,hLine4] = plotyy([blood.wave blood.wave],[blood.muhhb blood.muhbo2],[wl,wl],[cu_mvec(:,1),ni_mvec(:,1)],'semilogy');
set(hAx3,'xlim',[550 1000])


% figure
% [hAx,hLine1,hLine2] = plotyy(blood.wave,[blood.muhhb blood.muhbo2],wl,[cu_mvec(:,1),ni_mvec(:,1)]);
% set(hAx,'xlim',[500 1000],'ylim1',[0 2E5],'ylim2',[0 3])

[~,muap_hbo2,muap_hb,wlb]=myBlood(80,150,1);
%[~,~,muap_hb,wlb]=myBlood(0,150,1);

k_cu=muap_hbo2(wlb==800)/cu_mvec(wl==800,1);
k_ni=muap_hb(wlb==800)/ni_mvec(wl==800,1);


figure;
subplot(2,1,1)
plot(wlb,muap_hbo2,'r','linewidth',2)
hold on
plot(wlb,muap_hb,'color',[0.7 0 0.7]*1.2,'linewidth',2)
ylabel('\mu_a [mm^{-1}]')
set(gca,'xlim',[600 1000],'ytick',[0:0.2:0.6],'ylim',[0 0.6],'linewidth',2,'box','off','fontsize',fs)
l=legend('HbO_2, 150 g/L, SaO_2=80%','HHb')
set(l,'fontsize',fs-5,'linewidth',1)
grid on
subplot(2,1,2)
plot(wl,cu_mvec(:,1)*k_cu,'color',[0.4 0.4 1],'linewidth',2)
hold on
plot(wl,ni_mvec(:,1)*k_ni,'color',[0.2 0.5 0.2]*1.5,'linewidth',2)
xlabel('Wavelength [nm]')
ylabel('\mu_a [mm^{-1}]')
set(gca,'xlim',[600 1000],'ytick',[0:0.2:0.6],'ylim',[0 0.6],'linewidth',2,'box','off','fontsize',fs)
l=legend(['CuSO_4, ' num2str(round(k_cu,3)) 'M'],['NiSO_4, ' num2str(round(k_ni,3)) 'M'])
set(l,'fontsize',fs-5,'linewidth',1)
grid on