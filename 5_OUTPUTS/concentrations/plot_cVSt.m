clear all
clc

%INPUTS
pathcbtcDB = 'D:\Box\batch_new\Plots\'; %where to find the database with all cbtcs (generated with get_cbtc_db.m)
namecbtcDB = 'cbtcDB'; %base name of the cbtcs database

ana = 1; %1: analyze extraction rates; 2: screen depths
vx  = 0.0318/0.3; %upscaled vertical velocity

%OUTPUTS
ifig  = 3; 
nval  = 100;

tvec0 = linspace(0,400,400+1); %time vector
tvec  = linspace(0,400,nval+1); %time vector
cvec  = linspace(0,1,nval+1); %normalized concentration vector

%FOR RATES
if ana==1
    nrate = 4; irate0 = 1;
    ndpth = 1;
    nvar  = 4;   %number of analysis
    ncbtc = [147 150 150 150]; %number of cbtcs
    
    Lgd = ['Q_{out} = 6000.0 m^{3}/d'; 'Q_{out} = 3000.0 m^{3}/d'; 'Q_{out} = 1500.0 m^{3}/d'; 'Q_{out} = 750.0 m^{3}/d '];
    titlepos = [50 0.9 1.0];
    % Ind = ['a';'b';'c';'d'];

    p(1,:) = [0.12 0.53 0.35 0.4];
    p(2,:) = [0.50 0.53 0.35 0.4];
    p(3,:) = [0.12 0.10 0.35 0.4];
    p(4,:) = [0.50 0.10 0.35 0.4];

    cbPosition = [0.88 0.76 0.03 0.17];
    
    dca(1,:) = [5300.0 7349.40 9398.81];
    dca(2,:) = [5300.0 6395.33 7490.67];
    dca(3,:) = [5300.0 5884.34 6468.68];
    dca(4,:) = [5300.0 5633.75 5967.51];
    
elseif ana==2
    nrate = 2; irate0 = 2;
    ndpth = 3;
    nvar  = 3;   %number of analysis
    ivarv = [2 3 1];
    ncbtc = [150 150 150]; %number of cbtcs
    
    Lgd = ['d_{ts} = 50.0 m '; 'd_{ts} = 100.0 m'; 'd_{ts} = 150.0 m'];
    titlepos = [50 0.9 1.0];
    % Ind = ['a';'b';'c'];

    p(1,:) = [0.12 0.55 0.35 0.4];
    p(2,:) = [0.50 0.55 0.35 0.4];
    p(3,:) = [0.31 0.10 0.35 0.4];
    
    cbPosition = [0.69 0.30 0.03 0.17];
    
    dca(1,:) = [2650.00 3779.18 4908.36];
    dca(2,:) = [5300.00 6395.33 7490.67];
    dca(3,:) = [7950.00 9003.00 10056.00];

end

%create figure
if ishandle(ifig); close(ifig); end
hFig = figure(ifig); clf
set(gcf,'color','w');
set(hFig, 'Position', [100 30 1000 800])
Hprob = gobjects(nvar); %create series of subplots

%loop over rates and depths
for irate=irate0:nrate
    for idpth=1:ndpth
        if ana == 1
            ivar = irate;
        elseif ana == 2
            ivar = ivarv(idpth);
        end

        %ANALYSIS
        %get data
        filecbtcDB = [pathcbtcDB,namecbtcDB,'_d',int2str(idpth),'_r',int2str(irate),'.csv'];
        cbtcDB = csvread(filecbtcDB);

        %normalize concentrations
        cbtcDBnorm = zeros(length(tvec0),ncbtc(ivar));
        for i=1:ncbtc(ivar)
            cbtcDBnorm(:,i) = cbtcDB(:,i)/cbtcDB(length(tvec0),i);
        end

        %resample normalize concentrations (in time)
        cbtcDBnormS = zeros(length(tvec),ncbtc(ivar));
        for i=1:ncbtc(ivar)
            cbtcDBnormS(:,i) = interp1(tvec0,cbtcDBnorm(:,i),tvec);
        end

        %get probabilities to exceed concentrations
        cbtcprob = zeros(length(cvec),length(tvec));
        for i=1:length(tvec)
            for j=1:length(cvec)
                cbtcprob(j,i) = sum(cbtcDBnormS(i,:)>cvec(j))/ncbtc(ivar);
            end
        end

        %plot
        tg = zeros(length(tvec),length(tvec));
        for i=1:length(tvec);
            tg(:,i) = tvec(i);
        end
        cg = zeros(length(cvec),length(cvec));
        for i=1:length(cvec);
            cg(i,:) = cvec(i);
        end

        Hprob(ivar) = subplot(2,2,ivar);
        set(Hprob(ivar),'pos',p(ivar,:)); %.. subplot position
        h = pcolor(tg,cg,cbtcprob);
        colormap(parula)
        caxis([0 1])
        set(h, 'EdgeColor', 'none');
        %axis
        ax = gca;
        ax.FontSize = 16;
        ax.Layer = 'top';
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.YAxis.MinorTickValues = 0:0.1:1;
        ax.LineWidth = 1.0;
        ax.TickLength = [0.02 0.02];
        
        hold on
        [C,h] = contour(tg, cg, cbtcprob, [0.1 0.1], 'LineColor', [0.7 0.7 0.7], 'LineWidth', 1.5, 'LineStyle', '--', 'Showtext', 'on');
        clabel(C,h,'FontSize',10,'Color',[0.9 0.9 0.9])
        [C,h] = contour(tg, cg, cbtcprob, [0.5 0.5], 'LineColor', [0.7 0.7 0.7], 'LineWidth', 1.5, 'LineStyle', '--', 'Showtext', 'on');
        clabel(C,h,'FontSize',10,'Color',[0.5 0.5 0.5])
        [C,h] = contour(tg, cg, cbtcprob, [0.9 0.9], 'LineColor', [0.7 0.7 0.7], 'LineWidth', 1.5, 'LineStyle', '--', 'Showtext', 'on');
        clabel(C,h,'FontSize',10,'Color',[0.5 0.5 0.5])
        %hold off

        if ana==1; if ivar==1 || ivar==2; set(gca,'Xticklabel',[]); end; end %.. no x-axis label
        if ivar==2 || ivar==4; set(gca,'Yticklabel',[]); end %.. no x-axis label
        if ivar==1 || ivar==3; ylabel('c/c_0 [-]','FontSize',16); end
        if ivar==3 || ivar==4; xlabel('time [yr]','FontSize',16); end
        if ivar==2
            cb = colorbar; 
            cb.Position = cbPosition;
            cb.FontSize = 12;
        end
        text(220,0.1,Lgd(ivar,:),'Color',[0.5 0.5 0.5],'FontSize',12.7,'FontWeight','bold')
        
        %add characteristic times
        tca(ivar,:) = dca(ivar,:)/vx/365;
        plot1 = plot([tca(ivar,1) tca(ivar,3)],[0 1],'Color',[255/255 92/255 51/255],'LineStyle','--','LineWidth',4);
        plot1.Color(4) = 0.5;
    end
end

%print -dpng -r800 cdf_cVSt_rates_v4