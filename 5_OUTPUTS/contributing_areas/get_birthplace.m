%plot as function of distance from well
clear all

%WELL
nwll = 3;
%..#1
well(1,1) = 16560.0;
well(2,1) = 4550.0;
%..#2
well(1,2) = 16560.0;
well(2,2) = 3050.0;
%..#3
well(1,3) = 16560.0;
well(2,3) = 1550.0;

%Domain discretization
dx = 80.0;
dy = 50.0;
nx = 240;
ny = 120;

%2D (top surface) GRID
Lx = nx*dx;
Ly = ny*dy;
xg = dx/2:dx:Lx-dx/2; %center of cell
yg = dy/2:dy:Ly-dy/2; %center of cell
[Xg,Yg] = meshgrid(xg,yg);

%2D GRID / plot as distance
xgp = -(Lx-well(1,1))+dx/2 : dx : well(1,1)-dx/2; %center of cell
ygp = -(Ly/2)+dy/2 : dy : (Ly/2)-dy/2; %center of cell 
[Xgp,Ygp] = meshgrid(xgp,ygp);

%REAL
nreal = 1;
idpth = 1;
irate = 1;

tlimprob = [25.0 50.0 75.0 100.0 125.0 150.0 175.0 200.0 225.0 250.0 275.0 300.0 325.0 350.0];

%CBTC
pathCBTCs = 'D:\Box\batch_new\3_OUTPUTS_mnw\cbtc\homo\';
basename  = 'cbtc';

%BIRTHPLACE
pathBP = 'D:\Box\batch_new\3_OUTPUTS_mnw\part_id0\homo\';
basenameBP = 'rw3d_part_evol';

%preallocate variables
npartgBP = zeros(size(yg,2),size(xg,2));
npartg   = zeros(size(yg,2),size(xg,2));
npart_w1 = zeros(size(yg,2),size(xg,2));
npart_w2 = zeros(size(yg,2),size(xg,2));
npart_w3 = zeros(size(yg,2),size(xg,2));
tavpartg = zeros(size(yg,2),size(xg,2));
st_tg    = zeros(size(yg,2),size(xg,2));
cv_tg    = zeros(size(yg,2),size(xg,2));
kk       = zeros(size(yg,2),size(xg,2));
tpartg   = zeros(20000,size(yg,2),size(xg,2));
cumprob  = zeros(size(yg,2),size(xg,2),size(tlimprob,2));

%--------------------------------------------------------------------------
%ANALYSE DATA

for ireal=1:nreal
    %print realization on command window
    if ireal==1; fprintf('%s\n','working on realization:'); end
    fprintf('%d%s',ireal,'..');
    if ireal==nreal; fprintf('\n'); end
    
    if ireal==21 && idpth==1 && irate==1
        test = 0;
    else
    
    %get BIRTHPLACES
    fileBP = [pathBP,basenameBP,int2str(ireal),'_d',int2str(idpth),'_r',int2str(irate),'.dat'];
    BP     = dlmread(fileBP,' ', 3, 0);
    %..reorganize matrix
    sizej = size(BP,1);
    sizei = size(BP,2);
    BPb   = zeros(sizej,3);

    for j=1:sizej
        ii=1;
        for i=1:sizei
            if BP(j,i)~=0
                BPb(j,ii)=BP(j,i);
                %locate particle in grid 
                if ii==2
                    igrid = ceil(BPb(j,1)/dx);
                    jgrid = ceil(BPb(j,2)/dy);
                    npartgBP(jgrid,igrid) = npartgBP(jgrid,igrid)+1;
                end
                ii=ii+1;
            end

        end
    end
    clear BP
    
    %get CBTC
    filename = [pathCBTCs,basename,'_real',int2str(ireal),'_d',int2str(idpth),'_r',int2str(irate),'.dat'];
    fid      = fopen(filename);
    fgets(fid); fgets(fid);

    iwll = 0;
    kk   = zeros(size(yg,2),size(xg,2));
    while 1 %~feof(fid)
        data = fgets(fid);
        if ischar(data)
            testZone = sscanf(data, '%c', 4);
            if strcmp(testZone,'ZONE') == 1 %change TP zone => next well
                iwll = iwll+1;
            else
                dataspl = sscanf(data, '%f%f%d%*f%*f%*f', 131); %either data or blank line
                if size(dataspl) == [3 1]  %if data is [t,c,id]

                    %get location in grid and add arrival time
                    igrid = ceil( BPb(dataspl(3),1) /dx );
                    if iwll==1
                        jgrid = ceil( (BPb(dataspl(3),2)-1500.0) /dy );
%                         if jgrid>0 && jgrid<=ny
                            npart_w1(jgrid,igrid) = npart_w1(jgrid,igrid)+1;
%                         end
                    elseif iwll == 2
                        jgrid = ceil( (BPb(dataspl(3),2)) /dy );
%                         if jgrid>0 && jgrid<=ny 
                            npart_w2(jgrid,igrid) = npart_w2(jgrid,igrid)+1;
%                         end
                    elseif iwll == 3
                        jgrid = ceil( (BPb(dataspl(3),2)+1500.0) /dy );
%                         if jgrid>0 && jgrid<=ny 
                            npart_w3(jgrid,igrid) = npart_w3(jgrid,igrid)+1;
%                         end
                    end
%                     if jgrid>0 && jgrid<=ny 
                        kk(jgrid,igrid) = kk(jgrid,igrid)+1;
                        npartg(jgrid,igrid) = npartg(jgrid,igrid)+1;
                        tavpartg(jgrid,igrid) = tavpartg(jgrid,igrid) + dataspl(1);
                        tpartg(kk(jgrid,igrid),jgrid,igrid) = dataspl(1);
                        for it=1:size(tlimprob,2)
                            if (dataspl(1)/365.0)<tlimprob(it)
                                cumprob(jgrid,igrid,it) = cumprob(jgrid,igrid,it)+1;
                            end
                        end
%                     end
                end

            end %strcmp (next well)

        else
            break
        end %ischar

    end %while
    fclose(fid);
    end
end

%get probability to reach a well
prob = zeros(size(yg,2),size(xg,2));
for igrid=1:size(xg,2)
    for jgrid=1:size(yg,2)
        prob(jgrid,igrid) =  (npart_w1(jgrid,igrid)/npartgBP(jgrid,igrid)+npart_w2(jgrid,igrid)/npartgBP(jgrid,igrid)+npart_w3(jgrid,igrid)/npartgBP(jgrid,igrid))/nwll;
    end
end

%colormap
M = 72;
G = fliplr(linspace(0,1,M)) .';
myGmap = horzcat(G, G, G);

%get arrival time average and st. deviation per cell
for igrid=1:size(xg,2)
    for jgrid=1:size(yg,2)

        if (~isnan(npartg(jgrid,igrid)))
            %average
            tavpartg(jgrid,igrid)=tavpartg(jgrid,igrid)/npartg(jgrid,igrid)/365;
            
            %coefficent of variation
            sumtg = 0;
            for i=1:npartg(jgrid,igrid) 
                sumtg = sumtg+(( tpartg(i,jgrid,igrid)/365-tavpartg(jgrid,igrid) )^2);
            end
            st_tg(jgrid,igrid) = sqrt(1/npartg(jgrid,igrid)*sumtg);
            cv_tg(jgrid,igrid) = st_tg(jgrid,igrid)/tavpartg(jgrid,igrid);
            
            %cumulative prob arrival time
            for it=1:size(tlimprob,2)
                cumprob(jgrid,igrid,it) = cumprob(jgrid,igrid,it)/npartg(jgrid,igrid);
            end
        end

    end
end

%PLOTS

%reorganize grid -> distance from a well
prob2     = zeros(size(yg,2),size(xg,2));
tavpartg2 = zeros(size(yg,2),size(xg,2));
cv_tg2    = zeros(size(yg,2),size(xg,2));
cumprob2  = zeros(size(yg,2),size(xg,2),size(tlimprob,2));

for igrid=1:size(xg,2)
    for jgrid=1:size(yg,2)
        prob2(jgrid,igrid) = prob(size(yg,2)-jgrid+1,size(xg,2)-igrid+1);
        if isnan(prob2(jgrid,igrid)) 
            prob2(jgrid,igrid) = 0;
        end
        tavpartg2(jgrid,igrid) = tavpartg(size(yg,2)-jgrid+1,size(xg,2)-igrid+1);
        if isnan(tavpartg2(jgrid,igrid)) 
            tavpartg2(jgrid,igrid) = 0;
        end
        cv_tg2(jgrid,igrid) = cv_tg(size(yg,2)-jgrid+1,size(xg,2)-igrid+1);
        if isnan(cv_tg2(jgrid,igrid)) 
            cv_tg2(jgrid,igrid) = 0;
        end
        for it=1:size(tlimprob,2)
            cumprob2(jgrid,igrid,it) = cumprob(size(yg,2)-jgrid+1,size(xg,2)-igrid+1,it);
            if isnan(cumprob2(jgrid,igrid,it)) 
                cumprob2(jgrid,igrid,it) = 0;
            end
        end
    end
end



%--------------------------------------------------------------------------
%EXPORT TO FILE
pathfile = 'D:\Box\batch_new\Plots\Outputs\all_stats\';
filename = [pathfile,'stats','_d',int2str(idpth),'_r',int2str(irate),'.xlsx'];
xlswrite(filename,prob2,'prob')
xlswrite(filename,tavpartg2,'av_trav')
xlswrite(filename,cv_tg2,'cv_trav')
for it=1:size(tlimprob,2)
    xlswrite(filename,cumprob2(:,:,it),['cumprob_t',num2str(tlimprob(it))])
end

fclose('all');



%plot prob. of reaching well
figure(2*irate+nwll*(irate-1)+1)
clf
prob2(prob2 == 0) = NaN;
%npartg2 = npartg2./npartgBP; %prob. of reaching well 
h = pcolor(Xgp,Ygp,prob2(:,:));
caxis([0 0.1])
set(h, 'EdgeColor', 'none');
colormap(jet)   %colormap(myGmap)   %colormap(flipud(bone))

set(gca,'FontSize',16)
xlabel('distance from well in x-direction [m]','FontSize',16);
ylabel('distance from well in y-direction [m]','FontSize',16);
colorbar

% figout = ['Prob_reachwll_d',int2str(idpth),'_r',int2str(irate)];
% print(figout,'-djpeg')

hold on
%scatter(well(1,2),well(2,2),'x')


%plot avarage of arrival times per grid cell
figure(2*irate+nwll*(irate-1)+2)
clf
tavpartg2(tavpartg2 == 0) = NaN;
h = pcolor(Xgp,Ygp,tavpartg2);
caxis([0 400])
set(h, 'EdgeColor', 'none');
colormap(parula)

% strmin = 'Mean arrival times';
% text(-1000,2250,strmin,'FontSize',16);
set(gca,'FontSize',16)
xlabel('distance from well in x-direction [m]','FontSize',16);
ylabel('distance from well in y-direction [m]','FontSize',16);
colorbar

% figout = ['Mean_arrivals_d',int2str(idpth),'_r',int2str(irate)];
% print(figout,'-djpeg')

hold on
%scatter(well(1,2),well(2,2),'x')


% %plot cv of arrival times per grid cell
% figure(2*irate+nwll*(irate-1)+3)
% clf
% st_tg(st_tg == 0) = NaN;
% cv_tg2(cv_tg2 == 0) = NaN;
% %caxis([0 1.1])
% h = pcolor(Xgp,Ygp,cv_tg2);
% set(h, 'EdgeColor', 'none');
% colormap(parula)
% 
% % strmin = 'Coef. of Variation of arrival times';
% % text(-1000,2250,strmin,'FontSize',16);
% set(gca,'FontSize',16)
% xlabel('distance from well in x-direction [m]','FontSize',16);
% ylabel('distance from well in y-direction [m]','FontSize',16);
% colorbar
% 
% figout = ['CV_arrivals_d',int2str(idpth),'_r',int2str(irate)];
% print(figout,'-djpeg')
% 
% hold on
% %scatter(well(1,2),well(2,2),'x')

%plot cum. prob. arrival times
% for it=1:size(tlimprob,2)
%     figure(2*irate+nwll*(irate-1)+3+it)
%     clf
%     cumprob(cumprob == 0) = NaN;
%     h = pcolor(Xgp,Ygp,cumprob2(:,:,it));
%     caxis([0 1.0])
%     set(h, 'EdgeColor', 'none');
%     colormap(parula)
%     
%     strmin = ['Prob. arrival time < ',num2str(tlimprob(it)),' y'];
%     text(-1000,2250,strmin,'FontSize',16);
%     set(gca,'FontSize',16)
%     xlabel('distance from well in x-direction [m]','FontSize',16);
%     ylabel('distance from well in y-direction [m]','FontSize',16);
%     colorbar
%     
%     figout = ['Prob_r',int2str(irate),'2_t',num2str(tlimprob(it))];
%     print(figout,'-djpeg')
%     hold on
% end

