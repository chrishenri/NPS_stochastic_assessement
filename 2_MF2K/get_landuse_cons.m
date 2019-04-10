%generate land map
clear all
close all

nreal = 50;

%Output files
pathout_rchTR = 'D:\Transient_Recharge\0_MF2K_inputs\rch_TR_yearly_pack\';
pathout_rchSS = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\rch_SS_pack\';
pathout_rchSS_homo = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\rch_homo_pack\';
basenameout_rch = 'box_rch';

pathout_partSS = 'D:\Transient_Recharge\SteadyState\0_RW3D_inputs\npart_loc\';
pathout_partSS_homo = 'D:\Transient_Recharge\SteadyState\0_RW3D_inputs\npart_loc_homo\';
basenameout_part = 'npart_loc';

% pathout_cf = 'D:\Box\pesticide_proj\0_RW3D_inputs\cell_file\';
% basenameout_cf = 'cellfile_cons';
% 
% pathout_tf = 'D:\Box\pesticide_proj\0_RW3D_inputs\time_func\';
% basenameout_tf = 'tf_cons';

%2D domain
nx = 120;
ny = 60;
nz = 625;
dx = 160.0;
dy = 100.0;

%Hydrau. cond. matrix
pathKmat = 'D:\Box\batch_new\0_MF2K_inputs\K-mat\';
baseKmat = 'KmatMF';
K = [0.01 0.5 50 200]; %K values for mud, muddy sand, sand and gravel
%K = [200 50 0.5 0.01];  %K values for gravel, sand, muddy sand and mud

%Crops & fields
ncrops = 6;
crops  = ['Almond'; 'Citrus'; 'Corn  '; 'Cotton'; 'Grain '; 'Grapes'];
probcrp = [0.24 0.24 0.18 0.12 0.12 0.10];

field_ndx = 3;
field_ndy = 5;

ncat = 3;
cat = ['  Sand';'Medium';'  Clay'];

load('VZ_Output_NO3.mat')

%TIME
tdisc   = 3; %daily, monthly or yearly
tend    = 300; %in years

nvalt   = [7670 252 21]; %number of values in time-series (Hydrus simul)
yloop   = 10; %loop the n last year of the Hydrus time series until tend
itstart = [nvalt(tdisc)-yloop*365+1 nvalt(tdisc)-yloop*12+1 nvalt(tdisc)-yloop*1+1]; %when to start the loop

nvaltot = [tend*365+1 tend*12+1 tend+1];
sp_val  = linspace(0,tend*365,nvaltot(tdisc));


%Generate mesh
Lx = nx*dx;
Ly = ny*dy;
xg = dx/2:dx:Lx-dx/2; %center of cell
yg = dy/2:dy:Ly-dy/2; %center of cell
[Xg,Yg] = meshgrid(xg,yg);

nptot = 5E5;
%mtot  = 1E6;
mmult = 1E4;
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%REORGANIZE DATA

%clean sand, cumulative daily data
if tdisc == 1
    r_Crop_sand = zeros(nvalt(tdisc),ncrops);
    cons_Crop_sand = zeros(nvalt(tdisc),ncrops);
    j=1;
    for i=1:size(Almond.Sand.daily(:,1),1)
        t = rem(Almond.Sand.daily(i,1),j);
        if t==0
            r_Crop_sand(j,1) = Almond.Sand.daily(i,2);
            r_Crop_sand(j,2) = Citrus.Sand.daily(i,2);
            r_Crop_sand(j,3) = Corn.Sand.daily(i,2);
            r_Crop_sand(j,4) = Cotton.Sand.daily(i,2);
            r_Crop_sand(j,5) = Grain.Sand.daily(i,2);
            r_Crop_sand(j,6) = Grapes.Sand.daily(i,2);

            cons_Crop_sand(j,1) = Almond.Sand.daily(i,3);
            cons_Crop_sand(j,2) = Citrus.Sand.daily(i,3);
            cons_Crop_sand(j,3) = Corn.Sand.daily(i,3);
            cons_Crop_sand(j,4) = Cotton.Sand.daily(i,3);
            cons_Crop_sand(j,5) = Grain.Sand.daily(i,3);
            cons_Crop_sand(j,6) = Grapes.Sand.daily(i,3);

            j=j+1;
        end
    end
end

%GATHER DATA IN SINGLE VARIABLE
r_ = zeros(nvalt(tdisc),ncrops,3);
cons_ = zeros(nvalt(tdisc),ncrops,3);

for it=2:nvalt(tdisc)
    if tdisc == 1
        %RECHARGE RATES (m/d)
        %Sand
        r_(it,1,1) = (r_Crop_sand(it-1,1) - r_Crop_sand(it,1)) * 0.01;
        r_(it,2,1) = (r_Crop_sand(it-1,2) - r_Crop_sand(it,2)) * 0.01;
        r_(it,3,1) = (r_Crop_sand(it-1,3) - r_Crop_sand(it,3)) * 0.01;
        r_(it,4,1) = (r_Crop_sand(it-1,4) - r_Crop_sand(it,4)) * 0.01;
        r_(it,5,1) = (r_Crop_sand(it-1,5) - r_Crop_sand(it,5)) * 0.01;
        r_(it,6,1) = (r_Crop_sand(it-1,6) - r_Crop_sand(it,6)) * 0.01;
        %Medium
        r_(it,1,2) = (Almond.Med.daily(it-1,2)  - Almond.Med.daily(it,2)) * 0.01;
        r_(it,2,2) = (Citrus.Med.daily(it-1,2)  - Citrus.Med.daily(it,2)) * 0.01;
        r_(it,3,2) = (Corn.Med.daily(it-1,2)    - Corn.Med.daily(it,2))   * 0.01;
        r_(it,4,2) = (Cotton.Med.daily(it-1,2)  - Cotton.Med.daily(it,2)) * 0.01;
        r_(it,5,2) = (Grain.Med.daily(it-1,2)   - Grain.Med.daily(it,2))  * 0.01;
        r_(it,6,2) = (Grapes.Med.daily(it-1,2)  - Grapes.Med.daily(it,2)) * 0.01;
        %Clay
        r_(it,1,3) = (Almond.Clay.daily(it-1,2) - Almond.Clay.daily(it,2)) * 0.01;
        r_(it,2,3) = (Citrus.Clay.daily(it-1,2) - Citrus.Clay.daily(it,2)) * 0.01;
        r_(it,3,3) = (Corn.Clay.daily(it-1,2)   - Corn.Clay.daily(it,2))   * 0.01;
        r_(it,4,3) = (Cotton.Clay.daily(it-1,2) - Cotton.Clay.daily(it,2)) * 0.01;
        r_(it,5,3) = (Grain.Clay.daily(it-1,2)  - Grain.Clay.daily(it,2))  * 0.01;
        r_(it,6,3) = (Grapes.Clay.daily(it-1,2) - Grapes.Clay.daily(it,2)) * 0.01;

        %NO3 RATES (g/m2/d)
        %Sand
        cons_(it,1,1) = (cons_Crop_sand(it-1,1) - cons_Crop_sand(it,1)) / 0.0001 * 0.001;
        cons_(it,2,1) = (cons_Crop_sand(it-1,2) - cons_Crop_sand(it,2)) / 0.0001 * 0.001;
        cons_(it,3,1) = (cons_Crop_sand(it-1,3) - cons_Crop_sand(it,3)) / 0.0001 * 0.001;
        cons_(it,4,1) = (cons_Crop_sand(it-1,4) - cons_Crop_sand(it,4)) / 0.0001 * 0.001;
        cons_(it,5,1) = (cons_Crop_sand(it-1,5) - cons_Crop_sand(it,5)) / 0.0001 * 0.001;
        cons_(it,6,1) = (cons_Crop_sand(it-1,6) - cons_Crop_sand(it,6)) / 0.0001 * 0.001;
        %Medium
        cons_(it,1,2) = (Almond.Med.daily(it-1,3)  - Almond.Med.daily(it,3)) / 0.0001 * 0.001;
        cons_(it,2,2) = (Citrus.Med.daily(it-1,3)  - Citrus.Med.daily(it,3)) / 0.0001 * 0.001;
        cons_(it,3,2) = (Corn.Med.daily(it-1,3)    - Corn.Med.daily(it,3))   / 0.0001 * 0.001;
        cons_(it,4,2) = (Cotton.Med.daily(it-1,3)  - Cotton.Med.daily(it,3)) / 0.0001 * 0.001;
        cons_(it,5,2) = (Grain.Med.daily(it-1,3)   - Grain.Med.daily(it,3))  / 0.0001 * 0.001;
        cons_(it,6,2) = (Grapes.Med.daily(it-1,3)  - Grapes.Med.daily(it,3)) / 0.0001 * 0.001;
        %Clay
        cons_(it,1,3) = (Almond.Clay.daily(it-1,3) - Almond.Clay.daily(it,3)) / 0.0001 * 0.001;
        cons_(it,2,3) = (Citrus.Clay.daily(it-1,3) - Citrus.Clay.daily(it,3)) / 0.0001 * 0.001;
        cons_(it,3,3) = (Corn.Clay.daily(it-1,3)   - Corn.Clay.daily(it,3))   / 0.0001 * 0.001;
        cons_(it,4,3) = (Cotton.Clay.daily(it-1,3) - Cotton.Clay.daily(it,3)) / 0.0001 * 0.001;
        cons_(it,5,3) = (Grain.Clay.daily(it-1,3)  - Grain.Clay.daily(it,3))  / 0.0001 * 0.001;
        cons_(it,6,3) = (Grapes.Clay.daily(it-1,3) - Grapes.Clay.daily(it,3)) / 0.0001 * 0.001;

        
    elseif tdisc == 2
        %RECHARGE RATES (m/d)
        %Sand
        r_(it,1,1) = (Almond.Sand.monthly(it-1,2) - Almond.Sand.monthly(it,2)) * 0.01 / (365/12);
        r_(it,2,1) = (Citrus.Sand.monthly(it-1,2) - Citrus.Sand.monthly(it,2)) * 0.01 / (365/12);
        r_(it,3,1) = (Corn.Sand.monthly(it-1,2)   - Corn.Sand.monthly(it,2))   * 0.01 / (365/12);
        r_(it,4,1) = (Cotton.Sand.monthly(it-1,2) - Cotton.Sand.monthly(it,2)) * 0.01 / (365/12);
        r_(it,5,1) = (Grain.Sand.monthly(it-1,2)  - Grain.Sand.monthly(it,2))  * 0.01 / (365/12);
        r_(it,6,1) = (Grapes.Sand.monthly(it-1,2) - Grapes.Sand.monthly(it,2)) * 0.01 / (365/12);
        %Medium
        r_(it,1,2) = (Almond.Med.monthly(it-1,2)  - Almond.Med.monthly(it,2)) * 0.01 / (365/12);
        r_(it,2,2) = (Citrus.Med.monthly(it-1,2)  - Citrus.Med.monthly(it,2)) * 0.01 / (365/12);
        r_(it,3,2) = (Corn.Med.monthly(it-1,2)    - Corn.Med.monthly(it,2))   * 0.01 / (365/12);
        r_(it,4,2) = (Cotton.Med.monthly(it-1,2)  - Cotton.Med.monthly(it,2)) * 0.01 / (365/12);
        r_(it,5,2) = (Grain.Med.monthly(it-1,2)   - Grain.Med.monthly(it,2))  * 0.01 / (365/12);
        r_(it,6,2) = (Grapes.Med.monthly(it-1,2)  - Grapes.Med.monthly(it,2)) * 0.01 / (365/12);
        %Clay
        r_(it,1,3) = (Almond.Clay.monthly(it-1,2) - Almond.Clay.monthly(it,2)) * 0.01 / (365/12);
        r_(it,2,3) = (Citrus.Clay.monthly(it-1,2) - Citrus.Clay.monthly(it,2)) * 0.01 / (365/12);
        r_(it,3,3) = (Corn.Clay.monthly(it-1,2)   - Corn.Clay.monthly(it,2))   * 0.01 / (365/12);
        r_(it,4,3) = (Cotton.Clay.monthly(it-1,2) - Cotton.Clay.monthly(it,2)) * 0.01 / (365/12);
        r_(it,5,3) = (Grain.Clay.monthly(it-1,2)  - Grain.Clay.monthly(it,2))  * 0.01 / (365/12);
        r_(it,6,3) = (Grapes.Clay.monthly(it-1,2) - Grapes.Clay.monthly(it,2)) * 0.01 / (365/12);

        %SIMAZINE RATES (g/m2/d)
        %Sand
        cons_(it,1,1) = (Almond.Sand.monthly(it-1,3) - Almond.Sand.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,2,1) = (Citrus.Sand.monthly(it-1,3) - Citrus.Sand.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,3,1) = (Corn.Sand.monthly(it-1,3)   - Corn.Sand.monthly(it,3))   / 0.0001 * 0.001 / (365/12);
        cons_(it,4,1) = (Cotton.Sand.monthly(it-1,3) - Cotton.Sand.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,5,1) = (Grain.Sand.monthly(it-1,3)  - Grain.Sand.monthly(it,3))  / 0.0001 * 0.001 / (365/12);
        cons_(it,6,1) = (Grapes.Sand.monthly(it-1,3) - Grapes.Sand.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        %Medium
        cons_(it,1,2) = (Almond.Med.monthly(it-1,3)  - Almond.Med.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,2,2) = (Citrus.Med.monthly(it-1,3)  - Citrus.Med.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,3,2) = (Corn.Med.monthly(it-1,3)    - Corn.Med.monthly(it,3))   / 0.0001 * 0.001 / (365/12);
        cons_(it,4,2) = (Cotton.Med.monthly(it-1,3)  - Cotton.Med.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,5,2) = (Grain.Med.monthly(it-1,3)   - Grain.Med.monthly(it,3))  / 0.0001 * 0.001 / (365/12);
        cons_(it,6,2) = (Grapes.Med.monthly(it-1,3)  - Grapes.Med.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        %Clay
        cons_(it,1,3) = (Almond.Clay.monthly(it-1,3) - Almond.Clay.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,2,3) = (Citrus.Clay.monthly(it-1,3) - Citrus.Clay.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,3,3) = (Corn.Clay.monthly(it-1,3)   - Corn.Clay.monthly(it,3))   / 0.0001 * 0.001 / (365/12);
        cons_(it,4,3) = (Cotton.Clay.monthly(it-1,3) - Cotton.Clay.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        cons_(it,5,3) = (Grain.Clay.monthly(it-1,3)  - Grain.Clay.monthly(it,3))  / 0.0001 * 0.001 / (365/12);
        cons_(it,6,3) = (Grapes.Clay.monthly(it-1,3) - Grapes.Clay.monthly(it,3)) / 0.0001 * 0.001 / (365/12);
        
    elseif tdisc == 3
        %RECHARGE RATES (m/d)
        %Sand
        r_(it,1,1) = (Almond.Sand.yearly(it-1,2) - Almond.Sand.yearly(it,2)) * 0.01 / 365;
        r_(it,2,1) = (Citrus.Sand.yearly(it-1,2) - Citrus.Sand.yearly(it,2)) * 0.01 / 365;
        r_(it,3,1) = (Corn.Sand.yearly(it-1,2)   - Corn.Sand.yearly(it,2))   * 0.01 / 365;
        r_(it,4,1) = (Cotton.Sand.yearly(it-1,2) - Cotton.Sand.yearly(it,2)) * 0.01 / 365;
        r_(it,5,1) = (Grain.Sand.yearly(it-1,2)  - Grain.Sand.yearly(it,2))  * 0.01 / 365;
        r_(it,6,1) = (Grapes.Sand.yearly(it-1,2) - Grapes.Sand.yearly(it,2)) * 0.01 / 365;
        %Medium
        r_(it,1,2) = (Almond.Med.yearly(it-1,2)  - Almond.Med.yearly(it,2)) * 0.01 / 365;
        r_(it,2,2) = (Citrus.Med.yearly(it-1,2)  - Citrus.Med.yearly(it,2)) * 0.01 / 365;
        r_(it,3,2) = (Corn.Med.yearly(it-1,2)    - Corn.Med.yearly(it,2))   * 0.01 / 365;
        r_(it,4,2) = (Cotton.Med.yearly(it-1,2)  - Cotton.Med.yearly(it,2)) * 0.01 / 365;
        r_(it,5,2) = (Grain.Med.yearly(it-1,2)   - Grain.Med.yearly(it,2))  * 0.01 / 365;
        r_(it,6,2) = (Grapes.Med.yearly(it-1,2)  - Grapes.Med.yearly(it,2)) * 0.01 / 365;
        %Clay
        r_(it,1,3) = (Almond.Clay.yearly(it-1,2) - Almond.Clay.yearly(it,2)) * 0.01 / 365;
        r_(it,2,3) = (Citrus.Clay.yearly(it-1,2) - Citrus.Clay.yearly(it,2)) * 0.01 / 365;
        r_(it,3,3) = (Corn.Clay.yearly(it-1,2)   - Corn.Clay.yearly(it,2))   * 0.01 / 365;
        r_(it,4,3) = (Cotton.Clay.yearly(it-1,2) - Cotton.Clay.yearly(it,2)) * 0.01 / 365;
        r_(it,5,3) = (Grain.Clay.yearly(it-1,2)  - Grain.Clay.yearly(it,2))  * 0.01 / 365;
        r_(it,6,3) = (Grapes.Clay.yearly(it-1,2) - Grapes.Clay.yearly(it,2)) * 0.01 / 365;

        %SIMAZINE RATES (g/m2/d)
        %Sand
        cons_(it,1,1) = (Almond.Sand.yearly(it-1,3) - Almond.Sand.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,2,1) = (Citrus.Sand.yearly(it-1,3) - Citrus.Sand.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,3,1) = (Corn.Sand.yearly(it-1,3)   - Corn.Sand.yearly(it,3))   / 0.0001 * 0.001 / 365;
        cons_(it,4,1) = (Cotton.Sand.yearly(it-1,3) - Cotton.Sand.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,5,1) = (Grain.Sand.yearly(it-1,3)  - Grain.Sand.yearly(it,3))  / 0.0001 * 0.001 / 365;
        cons_(it,6,1) = (Grapes.Sand.yearly(it-1,3) - Grapes.Sand.yearly(it,3)) / 0.0001 * 0.001 / 365;
        %Medium
        cons_(it,1,2) = (Almond.Med.yearly(it-1,3)  - Almond.Med.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,2,2) = (Citrus.Med.yearly(it-1,3)  - Citrus.Med.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,3,2) = (Corn.Med.yearly(it-1,3)    - Corn.Med.yearly(it,3))   / 0.0001 * 0.001 / 365;
        cons_(it,4,2) = (Cotton.Med.yearly(it-1,3)  - Cotton.Med.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,5,2) = (Grain.Med.yearly(it-1,3)   - Grain.Med.yearly(it,3))  / 0.0001 * 0.001 / 365;
        cons_(it,6,2) = (Grapes.Med.yearly(it-1,3)  - Grapes.Med.yearly(it,3)) / 0.0001 * 0.001 / 365;
        %Clay
        cons_(it,1,3) = (Almond.Clay.yearly(it-1,3) - Almond.Clay.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,2,3) = (Citrus.Clay.yearly(it-1,3) - Citrus.Clay.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,3,3) = (Corn.Clay.yearly(it-1,3)   - Corn.Clay.yearly(it,3))   / 0.0001 * 0.001 / 365;
        cons_(it,4,3) = (Cotton.Clay.yearly(it-1,3) - Cotton.Clay.yearly(it,3)) / 0.0001 * 0.001 / 365;
        cons_(it,5,3) = (Grain.Clay.yearly(it-1,3)  - Grain.Clay.yearly(it,3))  / 0.0001 * 0.001 / 365;
        cons_(it,6,3) = (Grapes.Clay.yearly(it-1,3) - Grapes.Clay.yearly(it,3)) / 0.0001 * 0.001 / 365;

    end
end

tplot = linspace(1,nvalt(tdisc),nvalt(tdisc));
divt = [365 12 1];

hFigR = figure(1);
clf
set(gcf,'color','w');
set(hFigR, 'Position', [2500 50 800 800])
for icat=1:ncat
    hR(icat) = subplot(ncat,1,icat);
    for icrop=1:ncrops
        plot(tplot/divt(tdisc),r_(:,icrop,icat),'LineWidth',1.6)
        hold on
    end
    maxy = hR(icat).YLim(2);
    axis([0 nvalt(tdisc) 0 maxy])
    set(gca,'FontSize',16)
    xlabel('time [y]','FontSize',16);
    ylabel('recharge [m/m2/d]','FontSize',16);
    if icat==1
        legend(crops)
    end
end

hFigC = figure(2);
clf
set(gcf,'color','w');
set(hFigC, 'Position', [2500 50 800 800])
for icat=1:ncat
    hC(icat) = subplot(ncat,1,icat);
    for icrop=1:ncrops
        plot(tplot/divt(tdisc),cons_(:,icrop,icat),'LineWidth',1.6)
        hold on
    end
    maxy = hC(icat).YLim(2);
    axis([0 nvalt(tdisc) 0 maxy])
    set(gca,'FontSize',16)
    xlabel('time [y]','FontSize',16);
    ylabel('mass flux [g/m2/d]','FontSize',16);
    if icat==1
        legend(crops)
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%GET FIELDS MAP

% landuse = zeros(ny,nx,it);
% for it=1:nvaltot(tdisc) %FOR TRANSIENT FIELD SPATIAL DISTRIBUTION
    %randomly generate fields
    fields  = zeros(ny/field_ndy,nx/field_ndx);
    for ifield=1:nx/field_ndx
        for jfield=1:ny/field_ndy
            r = rand;
            fields(jfield,ifield) = sum(r >= cumsum([0, probcrp]));
        end
    end
    
    %distribute fields in grid
    landuse = zeros(ny,nx);
    for igrid=1:nx
        for jgrid=1:ny
            landuse(jgrid,igrid) = fields(ceil(jgrid/field_ndy),ceil(igrid/field_ndx));
        end
    end
    
    %plot fields
    figLand = figure(3);
    hR = pcolor(Xg,Yg,landuse(:,:));
    set(hR, 'EdgeColor', 'none');
    colormap(parula)
    colorbar
    set(gca,'FontSize',16)
    xlabel('x [m]','FontSize',16);
    ylabel('y [m]','FontSize',16);
    set(figLand,'Position', [2000, 500, 1200, 1200/3.2]);
    drawnow
    
    dlmwrite('D:\Transient_Recharge\landuse.dat',landuse)

% end

Header1_rch = '# MODFLOW2000 Recharge Package';
Header2_rch = 'PARAMETER  0';
Header3_rch = '         3        50';
Header4_rch = '         1         1';
Header5_rch = '        18   1.00000(10e12.4)                   -1     RECHARGE';

av_r = zeros(nvalt(tdisc),nreal);

mptot = zeros(nreal,1);
mftot = zeros(nreal,1);
%--------------------------------------------------------------------------
for ireal=1:nreal
    %----------------------------------------------------------------------
    %GET TRANSIENT MAPS
    nameMF = [pathKmat,baseKmat,'_',int2str(ireal),'.dat'];
    KmatMF   = dlmread(nameMF,'\t',[0 0 ny-1 nx-1]); %read only first layer (in MF format)
    Kmatb  = zeros(ny,nx);
    for iy=1:ny
        Kmatb(iy,:) = KmatMF(ny-iy+1,:);
    end
%     KmatMF(KmatMF==K(1)) = 1;
%     KmatMF(KmatMF==K(2)) = 2;
%     KmatMF(KmatMF==K(3)) = 2;
%     KmatMF(KmatMF==K(4)) = 3;
%     
%     Kmatb(Kmatb==K(1)) = 1;
%     Kmatb(Kmatb==K(2)) = 2;
%     Kmatb(Kmatb==K(3)) = 2;
%     Kmatb(Kmatb==K(4)) = 3;

    KmatMF(KmatMF==K(1)) = 3; %clay
    KmatMF(KmatMF==K(2)) = 2; %med
    KmatMF(KmatMF==K(3)) = 1; %sand
    KmatMF(KmatMF==K(4)) = 1; %sand
    
    Kmatb(Kmatb==K(1)) = 3;
    Kmatb(Kmatb==K(2)) = 2;
    Kmatb(Kmatb==K(3)) = 1;
    Kmatb(Kmatb==K(4)) = 1;
    
    
    %plot fields
%     fig = figure(6);
%     h = pcolor(Xg,Yg,Kmatb(:,:));
%     set(h, 'EdgeColor', 'none');
%     colormap(parula)
%     colorbar
%     set(gca,'FontSize',16)
%     xlabel('x [m]','FontSize',16);
%     ylabel('y [m]','FontSize',16);
%     set(fig,'Position', [2000, 500, 1200, 1200/3.2]);
%     drawnow
    
    %open file and generate rch package
    fileRCH_tr = [pathout_rchTR,basenameout_rch,'_',int2str(ireal),'.rch'];
    fileID_RCH_tr  = fopen(fileRCH_tr,'w');
    fprintf(fileID_RCH_tr,'%c',Header1_rch); fprintf(fileID_RCH_tr,'\n');
    fprintf(fileID_RCH_tr,'%c',Header2_rch); fprintf(fileID_RCH_tr,'\n');
    fprintf(fileID_RCH_tr,'%c',Header3_rch); fprintf(fileID_RCH_tr,'\n');
    
    %generate
    r_map    = zeros(ny,nx,nvalt(tdisc));
    cons_map = zeros(ny,nx,nvalt(tdisc));
    ncll_cons = 0;
    for it=itstart(tdisc):nvalt(tdisc)
        for igrid=1:nx
            for jgrid=1:ny
                r_map(jgrid,igrid,it) = r_(it,landuse(jgrid,igrid),KmatMF(jgrid,igrid));
                cons_map(jgrid,igrid,it) = cons_(it,landuse(jgrid,igrid),KmatMF(jgrid,igrid));
            end
        end
        av_r(it,ireal) = mean(mean(r_map(:,:,it)));
        if nnz(cons_map(:,:,it)) > ncll_cons
            ncll_cons = nnz(cons_map(:,:,it));
        end
        
        
        %PRINT RECHARGE MATRIX FOR STRESS PERIOD IT (MF2K RCH PACKAGE)
        fprintf(fileID_RCH_tr,'%c',Header4_rch); fprintf(fileID_RCH_tr,'\n');
        fprintf(fileID_RCH_tr,'%c',Header5_rch); fprintf(fileID_RCH_tr,'\n');
        %reshape recharge matrix (to vector)
        Rprint = reshape(r_map(:,:,it).', nx*ny, 1);
        ix=1;
        for i=1:nx*ny
            fprintf(fileID_RCH_tr,'%12.4e',Rprint(i,1));
            if ix==10
                fprintf(fileID_RCH_tr,'\n');
                ix=0;
            end
            ix = ix+1;
        end
        
        
%         for icrop=1:ncrops
%             for icat=1:ncat
%                 %PRINT FOR TR CONDITIONS
%                 %... GENERATE CELL FILE
%                 fileCF = [pathout_cf,basenameout_cf,'_crop',int2str(icrop),'_cat',int2str(icat),'_real',int2str(ireal),'.dat'];
%                 fileID_CF = fopen(fileCF,'w');
%                 fprintf(fileID_CF,'%d\n',ncll_cons);
%                 fprintf(fileID_CF,'%s\n','3');
%                 fprintf(fileID_CF,'%s\n','ix');
%                 fprintf(fileID_CF,'%s\n','iy');
%                 fprintf(fileID_CF,'%s\n','iz');
%                 
%                 %... WRITE CELLS
%                 for it=itstart(tdisc):nvalt(tdisc)
%                     for igrid=1:nx
%                     for jgrid=1:ny
%                         if cons_map(jgrid,igrid,it) > 0
%                             if print_conscll(jgrid,igrid) == 0
%                                 print_conscll(jgrid,igrid) = 1;
%                                 fprintf(fileID_CF,'%d\t%d\t%s\n',igrid,ny-jgrid+1,'1');
%                             end
%                         end
%                     end
%                     end
%                 end
%                 
%                 %... GENERATE TIME FUNCTION
%                 fileTF = [pathout_tf,basenameout_tf,'_crop',int2str(icrop),'_cat',int2str(icat),'_real',int2str(ireal),'.dat'];
%                 fileID_TF = fopen(fileTF,'w');
%                 fprintf(fileID_TF,'%s\n','Time function');
%                 fprintf(fileID_TF,'%s\n','inj_simazine');
%                 fprintf(fileID_TF,'%d\n',size(sp_val,2));
%                 
%                 %... WRITE TIME_FUNCTION
%                 it = itstart(tdisc);
%                 for its=1:size(sp_val,2)
%                     fprintf(fileID_TF,'%f\t%e\n',sp_val(its),cons_(it,icrop,icat)*dx*dy); 
%                     if it==nvalt(tdisc)
%                         it=itstart(tdisc);
%                     else
%                         it=it+1;
%                     end
%                 end
%             end
%         end

        
        
        %PLOT
%         figr = figure(2);
%         h = pcolor(Xg,Yg,r_map(:,:,it)); %wrong y-axis -> r_map is in MF format
%         set(h, 'EdgeColor', 'none');
%         %colormap(flipud(bone))
%         colormap(jet)
%         caxis( [ min(min(min(r_))) max(max(max(r_)))*0.7 ] )
%         colorbar
%         set(gca,'FontSize',16)
%         xlabel('x [m]','FontSize',16);
%         ylabel('y [m]','FontSize',16);
%         set(figr,'Position', [2000, 500, 1200, 1200/3.2]);
%         drawnow
%         figout = ['rechmap_t',int2str(it)];
%         print(figout,'-djpeg')
%         
%         
%         figsima = figure(3);
%         h = pcolor(Xg,Yg,sima_map(:,:,it)); %wrong y-axis -> r_map is in MF format
%         set(h, 'EdgeColor', 'none');
%         %colormap(flipud(bone))
%         colormap(jet)
%         caxis( [ min(min(min(sima_))) max(max(max(sima_)))*0.7 ] )
%         colorbar
%         set(gca,'FontSize',16)
%         xlabel('x [m]','FontSize',16);
%         ylabel('y [m]','FontSize',16);
%         set(figsima,'Position', [2000, 500, 1200, 1200/3.2]);
%         drawnow
%         figout = ['simamap_t',int2str(it)];
%         print(figout,'-djpeg')
%         
%         
%         figacet = figure(4);
%         h = pcolor(Xg,Yg,acet_map(:,:,it)); %wrong y-axis -> r_map is in MF format
%         set(h, 'EdgeColor', 'none');
%         %colormap(flipud(bone))
%         colormap(jet)
%         caxis( [ min(min(min(acet_))) max(max(max(acet_)))*0.7 ] )
%         colorbar
%         set(gca,'FontSize',16)
%         xlabel('x [m]','FontSize',16);
%         ylabel('y [m]','FontSize',16);
%         set(figacet,'Position', [2000, 500, 1200, 1200/3.2]);
%         drawnow
%         figout = ['acetmap_t',int2str(it)];
%         print(figout,'-djpeg')
    end
    fclose(fileID_RCH_tr);
    
    
    %open file and generate rch package for SS conditions
    fileRCH_ss = [pathout_rchSS,basenameout_rch,'_',int2str(ireal),'.rch'];
    fileID_RCH_ss  = fopen(fileRCH_ss,'w');
    fprintf(fileID_RCH_ss,'%c',Header1_rch); fprintf(fileID_RCH_ss,'\n');
    fprintf(fileID_RCH_ss,'%c',Header2_rch); fprintf(fileID_RCH_ss,'\n');
    fprintf(fileID_RCH_ss,'%c',Header3_rch); fprintf(fileID_RCH_ss,'\n');
    fprintf(fileID_RCH_ss,'%c',Header4_rch); fprintf(fileID_RCH_ss,'\n');
    fprintf(fileID_RCH_ss,'%c',Header5_rch); fprintf(fileID_RCH_ss,'\n');
    
    % get recharge rate averaged over the X considered years 
     
    sum_r = zeros(ny,nx);
    sum_cons = zeros(ny,nx);
    for it=itstart(tdisc):nvalt(tdisc)
        sum_r = sum_r + r_map(:,:,it);
        sum_cons = sum_cons + cons_map(:,:,it);
    end
    r_map_ss = sum_r / yloop;
    cons_map_ss = sum_cons / yloop;
    clear sum_r sum_cons
    
    % from MF format to Gslib format
    cons_map_ssGS = zeros(ny,nx);
    for iy=1:ny
        cons_map_ssGS(iy,:) = cons_map_ss(ny-iy+1,:);
    end
    cons_map_ss = cons_map_ssGS;
    clear cons_map_ssGS
    
    %reshape recharge matrix (to vector) and print
    Rprint = reshape(r_map_ss.', nx*ny, 1);
    ix=1;
    for i=1:nx*ny
        fprintf(fileID_RCH_ss,'%12.4e',Rprint(i,1));
        if ix==10
            fprintf(fileID_RCH_ss,'\n');
            ix=0;
        end
        ix = ix+1;
    end
    fclose(fileID_RCH_ss);
    
    %plot SS map
%     figRSS = figure(4);
%     set(figRSS,'Position', [2000, 500, 1200, 1200/3.2]);
%     hR = pcolor(Xg,Yg,r_map_ss(:,:));
%     set(hR, 'EdgeColor', 'none');
%     colormap(parula)
%     colorbar
%     set(gca,'FontSize',16)
%     xlabel('x [m]','FontSize',16);
%     ylabel('y [m]','FontSize',16);
    

    %GET RW3D INPUTS FOR PARTICLE INJECTIONS (CELL FILES AND TIME FUNCTIONS)
    %plot SS RECH map
%     figCSS = figure(5);
%     set(figCSS,'Position', [2000, 500, 1200, 1200/3.2]);
%     hR = pcolor(Xg,Yg,cons_map_ss(:,:));
%     set(hR, 'EdgeColor', 'none');
%     colormap(parula)
%     colorbar
%     set(gca,'FontSize',16)
%     xlabel('x [m]','FontSize',16);
%     ylabel('y [m]','FontSize',16);
    
    %get particle number map
    mftot(ireal,1) = sum(cons_map_ss(:));
    np_map = zeros(ny,nx);
    for igrid=1:nx
        for jgrid=1:ny
            np_map(jgrid,igrid) = floor(cons_map_ss(jgrid,igrid)*nptot/mftot(ireal,1));
        end
    end
%     figNPSS = figure(5);
%     set(figNPSS,'Position', [2000, 500, 1200, 1200/3.2]);
%     hR = pcolor(Xg,Yg,np_map(:,:));
%     set(hR, 'EdgeColor', 'none');
%     colormap(parula)
%     colorbar
%     set(gca,'FontSize',16)
%     xlabel('x [m]','FontSize',16);
%     ylabel('y [m]','FontSize',16);
    

    %open file and generate RW3D input for SS conditions
    nblocks = nx*ny;
    npinj   = sum(np_map(:));
    %mp      = mtot/npinj;
    mp_map  = zeros(ny,nx);
    filePART_ss = [pathout_partSS,basenameout_part,'_',int2str(ireal),'_ss','.dat'];
    fileID_PART_ss  = fopen(filePART_ss,'w');
    fprintf(fileID_PART_ss,'%d',nblocks);  fprintf(fileID_PART_ss,'\n');
    fprintf(fileID_PART_ss,'%d',npinj);    fprintf(fileID_PART_ss,'\n');
    %fprintf(fileID_PART_ss,'%10.6f',mp);   fprintf(fileID_PART_ss,'\n');
    fprintf(fileID_PART_ss,'%c','random'); fprintf(fileID_PART_ss,'\n');
    fprintf(fileID_PART_ss,'%c','random'); fprintf(fileID_PART_ss,'\n');
    fprintf(fileID_PART_ss,'%c','center'); fprintf(fileID_PART_ss,'\n');
    for igrid=1:nx
        for jgrid=1:ny
            if np_map(jgrid,igrid)>0
                mp_map(jgrid,igrid) = cons_map_ss(jgrid,igrid)/np_map(jgrid,igrid)*mmult;
            end
            fprintf(fileID_PART_ss,'%d\t%d\t%d\t%d\t%10.6f',igrid,jgrid,nz-1,np_map(jgrid,igrid),mp_map(jgrid,igrid));
            fprintf(fileID_PART_ss,'\n');
        end
    end
    
    mptot(ireal,1) = sum(mp_map(:));
    
    %HOMOGENEOUS CONDITIONS
    %open file and generate rch package for SS conditions
    fileRCH_ss_homo = [pathout_rchSS_homo,basenameout_rch,'_',int2str(ireal),'.rch'];
    fileID_RCH_ss_homo  = fopen(fileRCH_ss_homo,'w');
    fprintf(fileID_RCH_ss_homo,'%c',Header1_rch); fprintf(fileID_RCH_ss_homo,'\n');
    fprintf(fileID_RCH_ss_homo,'%c',Header2_rch); fprintf(fileID_RCH_ss_homo,'\n');
    fprintf(fileID_RCH_ss_homo,'%c',Header3_rch); fprintf(fileID_RCH_ss_homo,'\n');
    fprintf(fileID_RCH_ss_homo,'%c',Header4_rch); fprintf(fileID_RCH_ss_homo,'\n');
    fprintf(fileID_RCH_ss_homo,'%c',Header5_rch); fprintf(fileID_RCH_ss_homo,'\n');

    rch_av = mean(r_map_ss(:));
    ix=1;
    for i=1:nx*ny
        fprintf(fileID_RCH_ss_homo,'%12.4e',rch_av);
        if ix==10
            fprintf(fileID_RCH_ss_homo,'\n');
            ix=0;
        end
        ix = ix+1;
    end
    
    %open file and generate RW3D input for SS conditions  
    filePART_ss_homo = [pathout_partSS_homo,basenameout_part,'_',int2str(ireal),'_ss','.dat'];
    fileID_PART_ss_homo  = fopen(filePART_ss_homo,'w');
    
    av_mf = mean(cons_map_ss(:));
    av_np = floor(nptot/nblocks);
    mp_homo = av_mf/av_np*mmult;
    
    %np_av = floor(mean(np_map(:)));
    %npinj_homo = np_av*nblocks;
    %mp_homo = mptot(ireal,1)/npinj_homo;
    
    fprintf(fileID_PART_ss_homo,'%d',nblocks);      fprintf(fileID_PART_ss_homo,'\n');
    fprintf(fileID_PART_ss_homo,'%d',nptot);	fprintf(fileID_PART_ss_homo,'\n');
    %fprintf(fileID_PART_ss_homo,'%10.6f',mp_homo);	fprintf(fileID_PART_ss_homo,'\n');
    fprintf(fileID_PART_ss_homo,'%c','random');     fprintf(fileID_PART_ss_homo,'\n');
    fprintf(fileID_PART_ss_homo,'%c','random');     fprintf(fileID_PART_ss_homo,'\n');
    fprintf(fileID_PART_ss_homo,'%c','center');     fprintf(fileID_PART_ss_homo,'\n');
    
    for igrid=1:nx
        for jgrid=1:ny
            fprintf(fileID_PART_ss_homo,'%d\t%d\t%d\t%d\t%10.6f',igrid,jgrid,nz-1,av_np,mp_homo);
            fprintf(fileID_PART_ss_homo,'\n');
        end
    end
    
    clear r_map cons_map r_map_ss rch_av np_map np_av Rprint KmatMF Kmatb
    fclose('all');
end
