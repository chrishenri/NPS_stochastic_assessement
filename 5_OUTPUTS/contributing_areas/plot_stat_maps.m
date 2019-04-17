clear all

pathxls = 'D:\Box\batch_new\Plots\Outputs\all_stats\';
xlsnamebase = 'stats_v2';
sheetname_prob = 'prob';
sheetname_trav = 'av_trav';
sheetname_trcv = 'cv_trav';

pathwllinfo = 'D:\Box\batch_new\0_MF2K_inputs\mnw2_pack\';
namewllinfo = 'mc_scenarios';

ana = 2; %1 => analyze different rates; 2 => depths

qx = 0.0318;
qz = 0.0006;

if ana==1 %FOR RATES
    nvar  = 4; %number of analysis
    ir0   = 1; %start loop at irate=ir0
    nrate = 4;
    ndpth = 1;

    Lgd = ['Q_{out} = 6000.0 m^{3}/d'; 'Q_{out} = 3000.0 m^{3}/d'; 'Q_{out} = 1500.0 m^{3}/d'; 'Q_{out} = 750.0 m^{3}/d '];
    titlepos = [15000 -2700 1.0];
    Ind = ['a';'b';'c';'d'];

    p(1,:) = [0.15 0.75 0.8 0.22];
    p(2,:) = [0.15 (0.1+(0.75-0.1)/3*2) 0.8 0.22];
    p(3,:) = [0.15 (0.1+(0.75-0.1)/3) 0.8 0.22];
    p(4,:) = [0.15 0.1 0.8 0.22];

    cbpos(1,:) = [0.86 0.825 0.025 0.12];
    cbpos(2,:) = [0.86 (0.175+(0.825-0.175)/3*2) 0.025 0.12];
    cbpos(3,:) = [0.86 (0.175+(0.825-0.175)/3) 0.025 0.12];
    cbpos(4,:) = [0.86 0.175 0.025 0.12];
    
    dw = zeros(4,3);
    dw(:,1) = [100.0 100.0 100.0 100.0]; %depth of the top of the well screens

elseif ana==2 %FOR DEPTHS
    nvar  = 3; %number of analysis
    ir0   = 2; %start loop at irate=ir0
    nrate = 2;
    ndpth = 3;
    ivarv = [2 3 1];
    
    Lgd = ['d_{ts} = 50.0 m '; 'd_{ts} = 100.0 m'; 'd_{ts} = 150.0 m'];
    titlepos = [15000 -2700 1.0];
    Ind = ['a';'b';'c'];

    p(1,:) = [0.15 0.69 0.8 0.3];
    p(2,:) = [0.15 (0.1+(0.69-0.1)/2) 0.8 0.3];
    p(3,:) = [0.15 0.1 0.8 0.3];

    cbpos(1,:) = [0.86 0.80 0.025 0.15];
    cbpos(2,:) = [0.86 (0.21+(0.80-0.21)/2) 0.025 0.15];
    cbpos(3,:) = [0.86 0.21 0.025 0.15];
    
    dw = zeros(3,3);
    dw(:,1) = [50.0 100.0 150.0]; %depth of the top of the well screens

end

wellx = 16560.0; %location of extraction well (in x)

nQ = [625*2 313*2 157*2 79*2]; %number of cell needed for the total recharge to equal the extraction 

%close figures if open
ifig = [1 2];
if ishandle(ifig) 
    close(ifig)
end

%2D GRID / analysis
dx = 80.0;
dy = 50.0;
nx = 240;
ny = 120;
Lx = nx*dx;
Ly = ny*dy;

%2D GRID / plot as distance
xgp = -(Lx-wellx)+dx/2 : dx : wellx-dx/2; %center of cell
ygp = -(Ly/2)+dy/2 : dy : (Ly/2)-dy/2; %center of cell 
[Xgp,Ygp] = meshgrid(xgp,ygp);

Plim = zeros(nrate,ndpth);
%PLOT
hFig_prob = figure(ifig(1)); 
clf
set(gcf,'color','w');
set(hFig_prob, 'Position', [100 50 800 700])

hFig_trav = figure(ifig(2)); 
clf
set(gcf,'color','w');
set(hFig_trav, 'Position', [100 50 800 700])

dca = zeros(nvar,3);

Hprob = gobjects(nvar);
Htrav = gobjects(nvar);
for irate=ir0:nrate
    for idpth=1:ndpth
        if ana == 1
            ivar = irate;
        elseif ana == 2
            ivar = ivarv(idpth);
        end
        
        %Critical distance for no-pumping homogeneous case
        filewllinfo = [pathwllinfo,namewllinfo,'_d',int2str(idpth),'_r',int2str(irate),'.txt'];
        fid = fopen(filewllinfo,'rt');
        datawll=textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1,'CommentStyle','@');
        dw(ivar,2) = dw(ivar,1)+mean(datawll{1,7}(:))/2;
        dw(ivar,3) = dw(ivar,1)+mean(datawll{1,7}(:));
        fclose(fid);
        
        dca(ivar,:) = dw(ivar,:)/(qz/qx);
        
        
        %---------------------------
        %PLOT PROB
        hFig_prob = figure(ifig(1)); 
        
        %Get data
        prob = xlsread([pathxls,xlsnamebase,'_d',int2str(idpth),'_r',int2str(irate)],sheetname_prob);
        probV = sort(prob(:),'descend'); %sort all values of prob. matrix
        prob(prob == 0) = NaN;
        
        Hprob(ivar) = subplot(nvar,1,ivar);
        set(Hprob(ivar),'pos',p(ivar,:)); %.. subplot position
        hh = pcolor(Xgp,Ygp,prob(:,:));
        hold on
        %caxis([0 0.9])
        caxis([0 0.7792])
        set(hh, 'EdgeColor', 'none'); %.. no pcolor grid
        set(gca,'xgrid','on') %.. x axis grid
        set(gca,'FontSize',16) %.. axis font size
        if ivar~=nvar
            set(gca,'Xticklabel',[]) %.. no x-axis label
            %set(gca,'XColor',[1 0.1 0.1])
        end
        colormap(jet)
        
        %Legend
        cb = colorbar; 
        cb.Position = cbpos(ivar,:);
        cb.FontSize = 11;
        cb.Ticks = [0 0.25 0.5 0.75];
        %cb.Ticks = [0 0.45 0.9];
        
        %Get capture zone for an equivalent homogeneous case
        probHS = prob;
        Plim(irate,idpth) = probV(nQ(irate));
        probHS(probHS < Plim(irate,idpth)) = 0;
        probHS(probHS > 0) = Plim(irate,idpth);
        
        %Plot capture zone for an equivalent homogeneous case
        [CC,ii] = contourf(Xgp,Ygp,probHS(:,:)); %.. plot
        ii.Fill = 'off';
        ii.LineColor = [0.7 0.7 0.7];
        ii.LineWidth = 0.3;

        if ivar==nvar
            xlabel('distance from a well in the x-direction [m]','FontSize',15);
            t = text(-4700,700, 'distance from a well in the y-direction [m]','FontSize',15);
            set(t, 'rotation', 90)
        end
        %ylabel('distance from a well in the y-direction [m]','FontSize',16);
        hold on
        scatter(0,0,100,'x','MarkerEdgeColor',[0.8 0.0 0.0],'LineWidth',3)
        title(Lgd(ivar,:));
        set(get(gca,'title'),'Position',titlepos,'Color',[0.5 0.5 0.5],'FontSize',13.5)
        text(-1500,-2000,['(',Ind(ivar,:),')'],'Color',[0.5 0.5 0.5],'FontSize',14)
        xlim([-2000 18000])
        %xticks([0 4000 8000 12000 16000])
        
        %Plot critical distance for no-pumping homogeneous case
        scatter(dca(ivar,:),[0 0 0],75,'x','MarkerEdgeColor',[0.8 0.8 0.8],'LineWidth',2.5)
        plot(dca(ivar,:),[0 0 0],'Color',[0.8 0.8 0.8],'LineWidth',1.5)
        
        %---------------------------
        %PLOT TRAV
        hFig_trav = figure(ifig(2)); 
        
        %Get data average travel time
        trav = xlsread([pathxls,xlsnamebase,'_d',int2str(idpth),'_r',int2str(irate)],sheetname_trav);
        trav(trav == 0) = NaN;
        
        %Get data travel time coeff. of var. 
        trcv = xlsread([pathxls,xlsnamebase,'_d',int2str(idpth),'_r',int2str(irate)],sheetname_trcv);
        trcv(trcv == 0) = NaN;
        
        trsig = trcv.*trav;
        
        Htrav(ivar) = subplot(nvar,1,ivar);        
        set(Htrav(ivar),'pos',p(ivar,:)); %.. subplot position
        % hh = pcolor(Xgp,Ygp,trav(:,:));
        hh = pcolor(Xgp,Ygp,trcv(:,:));
        hold on
        % caxis([0 100])
        caxis([0 1.0])
        % caxis([0 400])
        set(hh, 'EdgeColor', 'none'); %.. no pcolor grid
        set(gca,'xgrid','on') %.. x axis grid
        set(gca,'FontSize',16) %.. axis font size
        if ivar~=nvar
            set(gca,'Xticklabel',[]) %.. no x-axis label
            %set(gca,'XColor',[1 0.1 0.1])
        end
        % colormap(flipud(jet))
        colormap(flipud(hot))
        
        %Legend
        cb = colorbar; 
        cb.Position = cbpos(ivar,:);
        cb.FontSize = 11;
        %cb.Ticks = [0 0.45 0.9];
        
        %plot hot spot
        % maxprob = max(prob(:));
        % probHS(probHS < 0.5*maxprob) = 0;
        % dataHS(dataHS < 0.25) = 0;
        
        %plot capture zone for an equivalent homogeneous case
        [CC,ii] = contourf(Xgp,Ygp,probHS(:,:)); %.. plot
        ii.Fill = 'off';
        ii.LineColor = [0.7 0.7 0.7];
        ii.LineWidth = 0.3;
        
        if ivar==nvar
            xlabel('distance from a well in the x-direction [m]','FontSize',15);
            t = text(-4700,700, 'distance from a well in the y-direction [m]','FontSize',15);
            %t = text(-4800,200, 'distance from a well in the y-direction [m]','FontSize',15);
            set(t, 'rotation', 90)
        end
        %ylabel('distance from a well in the y-direction [m]','FontSize',16);
        hold on
        scatter(0,0,100,'x','MarkerEdgeColor',[0.7 0.7 0.7],'LineWidth',3)
        title(Lgd(ivar,:));
        set(get(gca,'title'),'Position',titlepos,'Color',[0.5 0.5 0.5],'FontSize',13.5)
        text(-1500,-2000,['(',Ind(ivar,:),')'],'Color',[0.5 0.5 0.5],'FontSize',14)
        xlim([-2000 18000])
        %xticks([0 4000 8000 12000 16000])
        
        clear data dataHS
        
    end
end

% figure(ifig(1))
% print -depsc2 -painters prob_depths_v5
% print -dpng -r800 prob_rates_v5

% figure(ifig(2))
% print -depsc2 -painters av_trav_v4
% print -depsc2 -painters cv_trav_v4
% print -dpng -r800 av_trav_depths_v4
% print -dpng -r800 cv_trav_depths_v4

