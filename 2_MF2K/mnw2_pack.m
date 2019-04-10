%generate a series of mnw2 package, with different pumping rate, screen
%length, and top depth

%the well location is selected in order to always have 10ft of gravel/sand
%for each 100gpm of pumping rate. If this is not doable at the given well
%location, change this location until the criteria is fulfilled. 

clear all

% Fixed parameters
ndepth  = 3; %d in output file name
nrate   = 2; %r in output file name
nreal   = 50; %real in output file name

nt = 10;

pathout = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\mnw2_pack\';
outbase = 'box';
dx   = 160.0;
dy   = 100.0;
dz   = 0.4;
ncol = 120;
nrow = 60;
nlay = 625;

% Tsim file
path_ts  = 'D:\Box\batch_new\TPROGS\tsim\output\';
name_ts  = 'tsim_box.asc';
intGrav  = 1;
intSand  = 2;

K(1) = 200.0;
K(2) = 50.0;
K(3) = 0.5;
K(4) = 0.01;

% Initial well location. Can be modified if not enough sand/gravel at this
% location
nwll = 3; %number of well in a single realization
xwll = 16560.0;

wllname = ['W1';'W2';'W3'];

deltaxwll = 4*dx;
deltaywll = (dy*nrow)/(nwll+1);    
colmin = (xwll - deltaxwll/2)/dx;
colmax = (xwll + deltaxwll/2)/dx;
ywll   = zeros(nwll,1);
rowmin = zeros(nwll,1);
rowmax = zeros(nwll,1);
for iwll=1:nwll
    ywll(iwll,1)   = iwll*(dy*nrow)/(nwll+1);
    rowmin(iwll,1) = (ywll(iwll,1) - deltaywll/2)/dy;
    rowmax(iwll,1) = (ywll(iwll,1) + deltaywll/2)/dy;
end

nlaybel = 25; % number of layers bellow the bottom of the well
maxiter = 40; % max. number of iteration to find well location

% Pumping well parameters - could define distribution
depth_topwll(1) = 100.0;   % m
depth_topwll(2) = 150.0;   % m
depth_topwll(3) = 50.0;    % m

rate(1) = -6000.0;  % m3/d
rate(2) = -3000.0;  % m3/d
rate(3) = -1500.0;  % m3/d
rate(4) = -750.0; % m3/d

% File to keep track of all scenarios
basenametr  = 'mc_scenarios';

%-------------------------------------------------------------------------
% Build well package
for ireal=1:nreal
    if ireal==1; fprintf('%s\n','working on realization:'); end
    fprintf('%d%s',ireal,'..');
    if ireal==nreal; fprintf('\n'); end
    
    % .. get tsim file for this realization
    fname_ts = [path_ts,name_ts,int2str(ireal)];
    %ts = dlmread(fname_ts);
    fid = fopen(fname_ts,'rt');
    ts = textscan(fid,'%f','Headerlines',1);
    fclose(fid);
    
    tsmat = zeros(ncol,nrow,nlay);
    % reorganize ts
    i=1;
    for iz=1:nlay
        for iy=1:nrow
            for ix=1:ncol
                tsmat(ix,iy,iz) = ts{1,1}(i,1);
                i=i+1;
            end
        end
    end
    
    for irate=2:nrate
        thick_coarse_tot = 3.048*abs(rate(irate))/545.10;
        % get col,row
        col = zeros(nwll,1);
        row = zeros(nwll,1);        
        for iwll=1:nwll
            col(iwll,1) = int8(xwll/dx); % for intial well location
            row(iwll,1) = int8(ywll(iwll,1)/dy); % for intial well location
        end
        
        for idepth=1:ndepth
            fnametr = [basenametr,'_d',int2str(idepth),'_r',int2str(irate),'.txt'];
            if ireal==1
                fileIDtr = fopen(fnametr,'w');
                fprintf(fileIDtr,'     ireal      rate     depth      iwll       col       row    wll_th    coa_th \n');
            else
                fileIDtr = fopen(fnametr,'a');
            end
            
            % get laytopwel
            laytopwel = depth_topwll(idepth)/dz+1;
            
            % get laybotwel
            thick_coarse = zeros(nwll,1);
            sumK = zeros(nwll,1);
            print = zeros(nwll,1);
            dim_zwll = zeros(nwll,1);
            iter = zeros(nwll,1);
            
            for iwll=1:nwll
                converge = 0;
                
                while converge == 0
                    thick_coarse(iwll,1) = 0.0;
                    sumK(iwll,1) = 0.0;
                    for laybotwel = laytopwel:(nlay-nlaybel)
                        sumK(iwll,1) = sumK(iwll,1) + K(tsmat(col(iwll,1),row(iwll,1),laybotwel));
                        if tsmat(col(iwll,1),row(iwll,1),laybotwel) == 1 || tsmat(col(iwll,1),row(iwll,1),laybotwel) == 2
                            thick_coarse(iwll,1) = thick_coarse(iwll,1) + dz;
                        end

                        % if accumalated enough sand and gravel
                        if thick_coarse(iwll,1) >= thick_coarse_tot
                            print(iwll,1)  = 1;

                            % get layers with well
                            dim_zwll(iwll,1) = laybotwel-laytopwel+1;
                            %lay = zeros(dim_zwll(iwll,1),1);
                            for izwll=1:dim_zwll(iwll,1)
                                lay(izwll,iwll) = laytopwel+izwll-1;
                            end
                            converge = 1;
                            break
                        end
                    end
                    
                    if converge == 0    %didn't reach the required amount of coarse material -> change location
                        if iter(iwll,1) == maxiter
                            fprintf(fileIDtr,'%c','Could not find correct well location for ireal=');fprintf(fileIDtr,'%d',ireal); fprintf(fileID,'\n');
                            converge = 1;
                            %error('cannot find appropriate well location');
                        end

                        col(iwll,1) = round((colmax-colmin)*rand(1) + colmin); 
                        row(iwll,1) = round((rowmax(iwll,1)-rowmin(iwll,1))*rand(1) + rowmin(iwll,1));
                        iter(iwll,1) = iter(iwll,1) + 1;
                    end
                    
                end
                
            end
            
            
            %check if all wells are designed
            printf = 1;
            nwlldata = zeros(nwll,1);
            for iwll=1:nwll
                if print(iwll,1) == 0
                    printf = 0;
                end
            end
            
            % Write output (.wel package)
            if printf == 1
                % write well package
                Header1 = '        50   0';
                Header2 = 'THIEM 0 0 0 0';
                Header3 = '  1.000000e+000';
                Header4 = '                                Stress Period';
                var = 0;

                fname   = [pathout,outbase,'_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'.mnw2'];
                fileID  = fopen(fname,'w'); 
                fprintf(fileID,'%10d',nwll); fprintf(fileID,'%c',Header1); fprintf(fileID,'\n');
                
                nwlldata = zeros(nwll,1);
                for iwll=1:nwll
                    nwlldata(iwll,1) = nnz(lay(:,iwll));
                    fprintf(fileID,'%20s',wllname(iwll,:));fprintf(fileID,'%4d',nwlldata(iwll,1)); fprintf(fileID,'\n');
                    fprintf(fileID,'%c',Header2); fprintf(fileID,'\n');
                    fprintf(fileID,'%c',Header3); fprintf(fileID,'\n');
                    for izwll=1:dim_zwll(iwll,1)
                        fprintf(fileID,'%d\t%d\t%d\n',lay(izwll,iwll),row(iwll,1),col(iwll,1));
                    end
                                        
                    % Recap parameters in text file
                    fprintf(fileIDtr,'%10d%10.2f%10.3f%10d%10d%10d%10.3f%10.3f\n',ireal,rate(irate),depth_topwll(idepth),iwll,col(iwll,1),row(iwll,1),dim_zwll(iwll,1)*dz,thick_coarse(iwll,1));
                end
                for it=1:nt
                    numnum = numel(num2str(it));
                    fmt = ['%',int2str(numnum)+1,'d'];
                    fprintf(fileID,'%2d',nwll); fprintf(fileID,'%c',Header4); fprintf(fileID,fmt,it); fprintf(fileID,'\n');
                    for iwll=1:nwll
                        fprintf(fileID,'%20s',wllname(iwll,:)); fprintf(fileID,'%15e',rate(irate)); fprintf(fileID,'\n');
                    end
                end
                
                fclose(fileID);
            end
            
            clear thick_coarse sumK print dim_zwll lay ratelay
        end %depth
    end %rate

end %real
fclose('all');
