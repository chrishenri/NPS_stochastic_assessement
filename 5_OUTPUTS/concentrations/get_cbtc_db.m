clear all
clc

%INPUTS
pathCBTC = 'D:\Box\batch_new\3_OUTPUTS_mnw\cbtc\';
baseCBTC = 'cbtc';

%OUTPUTS
printout = 1; %1 if print the concentration database
pathfile = 'D:\Box\batch_new\Plots\';
filenameOUT = 'cbtcDB_d3_r2.csv';

nreal = 50;
realfail = [ ];

idpth = 3;
irate = 2;

tvec  = linspace(0,400,400+1); %time vector at which concentration are analyzed

%to convert cbtcs to concentration signals
Qout = [6000.0 3000.0 1500.0 750.0]; %extraction rate in m3/d
mmult = 1E4;
CellArea = 160*100;
mf_to_c = CellArea/mmult/Qout(irate);


dataB = zeros(length(tvec),150);
k = 0;
for ireal=1:nreal
    %print realization on command window
    if ireal==1; fprintf('%s','working on realization: '); end
    if ireal>1
    	for j=0:log10(ireal-1)
        	fprintf('\b'); % delete previous counter display
    	end
    end
    fprintf('%d', ireal);
    pause(.05); % allows time for display to update
    if ireal==nreal; fprintf('\n'); end
    
    if not(ismember(ireal,realfail))
    
     %-------------------------------------
    %CBTCs 
    filename = [pathCBTC,baseCBTC,'_real',int2str(ireal),'_d',int2str(idpth),'_r',int2str(irate),'.dat'];
    fid      = fopen(filename);
    fgets(fid); fgets(fid);
    
    %get cbtcs dimensions
    nbtc  = 0;
    while 1
        data = fgets(fid);
        if ischar(data)
            testZone = sscanf(data, '%c', 4);
            if strcmp(testZone,'ZONE') == 1 %change TP zone => next well
                nbtc = nbtc+1;
                ndata(nbtc) = 0;
            else
                ndata(nbtc) = ndata(nbtc)+1;
            end
        else
            break
        end
    end
    fclose('all');
    ndata = ndata-1;
    
    %get and analyze btcs
    for ibtc=1:nbtc
        if ibtc==1
            r1 = 1+2*ibtc;
            r2 = r1+ndata(ibtc)-1;
        else
            r1 = sum(ndata(1:ibtc-1))+1+2*ibtc;
            r2 = r1+ndata(ibtc)-1;
        end
        databtc = dlmread(filename,' ',[r1 1 r2 18]); 
        data2 = databtc(databtc ~= 0);
        clear databtc
        testdata = size(data2); testdata = testdata(1,1);
        if testdata>0
            data3(:,1) = data2(1:ndata(ibtc))./365; %t
            data3(:,2) = (data2(ndata(ibtc)+1:ndata(ibtc)*2)).*mf_to_c; %c in mg/L
        end
        clear data2
        
        %resample and gather data
        k = k+1;
        
        [timef, index] = unique(data3(:,1)); %get rid of non-unique points in times, if any
        concf = data3(index,2);
        dataB(:,k) = interp1(timef,concf,tvec);
        
        %dataB(:,k) = interp1(data3(:,1),data3(:,2),tvec);
        dataB(isnan(dataB))=0;
        [M,I] = max(dataB(:,k));
        dataB(I:length(tvec),k) = M;
        
        clear data3
%         figure(1)
%         plot(data3(:,1),data3(:,2))
%         hold on
%         plot(tvec,dataB(:,k),'LineStyle','--')
        
    end
    
    end
end
%export DB
if printout==1
csvwrite([pathfile,filenameOUT],dataB)
end
    