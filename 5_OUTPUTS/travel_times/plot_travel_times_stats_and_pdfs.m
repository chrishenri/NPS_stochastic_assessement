clear all

ifig = 4;
ana = 1; %1 for analyzing rates; 2 for depths

%---------------------------------------------
%---------------------------------------------
%PDFs
%rates analysis
% pathxls = 'D:\Box\batch_new\Plots\Outputs\pdfs\rate\';
% sheets = ['r750 '; 'r1500'; 'r3000'; 'r6000'];
% Lgd = ['Q_{out} = 750.0 m^{3}/d '; 'Q_{out} = 1500.0 m^{3}/d'; 'Q_{out} = 3000.0 m^{3}/d'; 'Q_{out} = 6000.0 m^{3}/d'];
% xlsname = 'all_rates_v2.xlsx';
% legpos = [0.3 0.74 0.2 0.1];
% 
% teq(4,:) = [136.99 189.96 242.93];
% teq(3,:) = [136.99 165.30 193.61];
% teq(2,:) = [136.99 152.09 167.19];
% teq(1,:) = [136.99 145.61 154.24];


pathxls = 'D:\Box\batch_new\Plots\Outputs\pdfs\dpth\';
sheets = ['d50 '; 'd100'; 'd150'];
Lgd = ['d_{ts} = 50.0 m '; 'd_{ts} = 100.0 m'; 'd_{ts} = 150.0 m'];
xlsname = 'all_depth_v2.xlsx';
ifig = 3;
legpos = [0.38 0.74 0.1 0.1];

teq(1,:) = [68.49	97.67	126.86];
teq(2,:) = [136.98	165.29	193.60];
teq(3,:) = [205.47	232.69	259.91];


cellsheets = cellstr(sheets);

for ish=1:size(cellsheets,1)
    sheetname = cellsheets{ish,1};
    data(:,:,ish) = xlsread([pathxls,xlsname],sheetname);
end

hFig = figure(ifig); 
clf
set(gcf,'color','w');
set(hFig, 'Position', [100 0 1200 800])

colors(:,1) = [255/255 212/255 58/255];
colors(:,2) = [174/255 218/255 226/255];
colors(:,3) = [237/255 115/255 111/255];
colors(:,4) = [98/255 188/255 146/255];


%--------------------------
%plot t5
H = subplot(3,2,1);
set(H,'pos',[0.1 0.66 0.45 0.25]); %.. subplot position
for ish=1:size(cellsheets,1)
    %plot pdf
    plot(data(:,2,ish),data(:,11,ish),'LineWidth',2,'Color',colors(:,ish))
    hold on
end

set(gca,'FontSize',16)
set(gca,'xgrid','on')
set(gca,'Xticklabel',[]) %.. no x-axis label

%     xlabel('arrival time [y]','FontSize',16);
ylabel('travel time pdf [1/yr]','FontSize',16);
% maxy = 0.015;
maxy = 0.02;
%axis([0 500 0 maxy])
axis([0 400 0 maxy])

title('t_{5}','FontSize',18)
% set(get(gca,'title'),'Position',[460 0.7*maxy 1.0],'Color',[0.5 0.5 0.5])
% text(450,0.15*maxy,'(a)','Color',[0.5 0.5 0.5],'FontSize',16)
set(get(gca,'title'),'Position',[370 0.6*maxy 1.0],'Color',[0.5 0.5 0.5])
text(15,0.8*maxy,'a)','Color',[0.5 0.5 0.5],'FontSize',16)
hold on

legend({Lgd},'FontSize',14,'Color',[1 1 1],'EdgeColor',[1 1 1],'position',legpos);    

%plot teq
for ish=1:size(cellsheets,1)
    plot1 = plot([teq(ish,1) teq(ish,1)],[0 maxy],'LineWidth',3,'LineStyle','--','Color',colors(:,ish));
    plot1.Color(4) = 0.3;
end

%--------------------------
%plot t50
H = subplot(3,2,3);
set(H,'pos',[0.1 0.38 0.45 0.25]); %.. subplot position
for ish=1:size(cellsheets,1)
    plot(data(:,5,ish),data(:,14,ish),'LineWidth',2,'Color',colors(:,ish))
    hold on
end

set(gca,'FontSize',16)
set(gca,'xgrid','on')
set(gca,'Xticklabel',[]) %.. no x-axis label

%     xlabel('arrival time [y]','FontSize',16);
ylabel('travel time pdf [1/yr]','FontSize',16);
maxy = 0.015;
%axis([0 500 0 maxy])
axis([0 400 0 maxy])

title('t_{50}','FontSize',18)
% set(get(gca,'title'),'Position',[460 0.7*maxy 1.0],'Color',[0.5 0.5 0.5])
% text(450,0.15*maxy,'(a)','Color',[0.5 0.5 0.5],'FontSize',16)
set(get(gca,'title'),'Position',[370 0.6*maxy 1.0],'Color',[0.5 0.5 0.5])
text(15,0.8*maxy,'b)','Color',[0.5 0.5 0.5],'FontSize',16)
hold on

%plot teq
for ish=1:size(cellsheets,1)
    plot1 = plot([teq(ish,2) teq(ish,2)],[0 maxy],'LineWidth',3,'LineStyle','--','Color',colors(:,ish));
    plot1.Color(4) = 0.3;
end

%--------------------------
%plot t90
H = subplot(3,2,5);
set(H,'pos',[0.1 0.1 0.45 0.25]); %.. subplot position
for ish=1:size(cellsheets,1)
    plot(data(:,8,ish),data(:,17,ish),'LineWidth',2,'Color',colors(:,ish))
    hold on
end
set(gca,'FontSize',16)
set(gca,'xgrid','on')

xlabel('travel time [yr]','FontSize',16);
ylabel('travel time pdf [1/yr]','FontSize',16);
maxy = 0.01;
%axis([0 500 0 maxy])
axis([0 400 0 maxy])

title('t_{95}','FontSize',18)
% set(get(gca,'title'),'Position',[460 0.7*maxy 1.0],'Color',[0.5 0.5 0.5])
% text(450,0.15*maxy,'(a)','Color',[0.5 0.5 0.5],'FontSize',16)
set(get(gca,'title'),'Position',[370 0.6*maxy 1.0],'Color',[0.5 0.5 0.5])
text(15,0.8*maxy,'c)','Color',[0.5 0.5 0.5],'FontSize',16)
hold on
% end

%plot teq
for ish=1:size(cellsheets,1)
    plot1 = plot([teq(ish,3) teq(ish,3)],[0 maxy],'LineWidth',3,'LineStyle','--','Color',colors(:,ish));
    plot1.Color(4) = 0.3;
end

clear all

%---------------------------------------------
%---------------------------------------------
%STATS
pathxls = 'D:\Box\batch_new\Plots\';
xlsname = 'stat_arrtime_v2.xlsx';

ana = 2;

sheets{1,:} = ['plot_mean_qout'; 'plot_cv_qout  ']; %rates
sheets{2,:} = ['plot_mean_dts'; 'plot_cv_dts  ']; %depths
Lgd = ['t_{5} '; 't_{50}'; 't_{95}'];


%position of the legend
legpos(1,:) = [0.72 0.435 0.1 0.15]; %rates
legpos(2,:) = [0.72 0.435 0.1 0.15]; %depths

%position of the title (a or b)
titlepos1(1,:) = [750 0.06*(300-50)+50 1.0]; %rates
titlepos1(2,:) = [140 0.07*400 1.0]; %depths

titlepos2(1,:) = [750 0.06*(0.4-0.2)+0.20 1.0];
titlepos2(2,:) = [140 0.13*(0.5-0.1)+0.1 1.0];

%plot label, x axis
xlab(1,:) = 'Q_{out} [m';%^{3} d^{-1}]';
xlab(2,:) = 'd_{ts} [m]';

%plot colors
colors(:,1) = [255/255 212/255 58/255];
colors(:,2) = [174/255 218/255 226/255];
colors(:,3) = [237/255 115/255 111/255];
colors(:,4) = [98/255 188/255 146/255];

%read data
cellsheets = cellstr(sheets{ana,:});

for ish=1:size(cellsheets,1)
    sheetname = cellsheets{ish,1};
    data(:,:,ish) = xlsread([pathxls,xlsname],sheetname);
end

%plot av(ti)
H = subplot(3,2,2);
set(H,'pos',[0.65 0.55 0.28 0.35]); %.. subplot position
for it=2:4
    plot(data(:,1,1),data(:,it,1),'-o','LineWidth',2,'MarkerEdgeColor','auto','MarkerFaceColor','auto','MarkerSize',7)
    hold on
end
set(gca,'FontSize',16)
set(gca,'ygrid','on')
set(gca,'Xticklabel',[]) %.. no x-axis label
ylabel('\langlet_i\rangle [yr]','FontSize',16);
title('d)','FontSize',16)
set(get(gca,'title'),'Position',titlepos1(ana,:),'Color',[0.5 0.5 0.5])
hold on

%plot std(ti)
H = subplot(3,2,6);
set(H,'pos',[0.65 0.12 0.28 0.35]); %.. subplot position
for it=2:4
    %plot(data(:,1,2),sqrt(data(:,it,2)),'-o','LineWidth',2,'MarkerFaceColor','auto','MarkerSize',7)
    plot(data(:,1,2),data(:,it,2),'-o','LineWidth',2,'MarkerFaceColor','auto','MarkerSize',7)
    hold on
end
set(gca,'FontSize',16)
set(gca,'ygrid','on')

xlabel(xlab(ana,:),'FontSize',16);
ylabel('CV_{t_{i}} [-]','FontSize',16);
title('e)','FontSize',16)
set(get(gca,'title'),'Position',titlepos2(ana,:),'Color',[0.5 0.5 0.5])
hold on

legend({Lgd},'FontSize',14,'Color',[1 1 1],'EdgeColor',[1 1 1],'position',legpos(ana,:),'Orientation','horizontal');

print -depsc2 -painters allstat_depths


