%plot nreal CBTCs and analyse results (pdf of t15 and t90)
clear all

pathCBTCs = 'D:\Box\batch_new\3_OUTPUTS_mnw\cbtc\';
basename  = 'cbtc';

pathwllinfo = 'D:\Box\batch_new\0_MF2K_inputs\mnw2_pack\';
wllinfonam  = 'mc_scenarios';

% pathCBTCs = 'D:\Box\batch_new\BW_simulation\3_OUTPUTS_fwadv\cbtc\';
% basename  = 'cbtc';

nreal = 50;
idpth = 1;
irate = 2;

nwell = 1;
iwllp = 2;

tana = [1 5 10 25 50 75 90 95];
win_min_ind = [2 3 4];
win_max_ind = [8 7 6];
ti   = zeros((nreal-1)*nwell,size(tana,2));
aw   = zeros((nreal-1)*nwell,size(win_min_ind,2));

tic
it=0;
for ireal=1:nreal
    if ireal==1
        fprintf('%s\n','working on realization:')
    end
    fprintf('%d%s',ireal,'..');
    if ireal==nreal
        fprintf('\n')
    end
    
%     if ireal==5||ireal==21||ireal==29||ireal==30||ireal==32||ireal==50   %realization didn't worked out irate 1
%     if ireal==50||ireal==27              %realization didn't worked out irate 2
    if ireal==50                         %realization didn't worked out irate 3
%     if ireal==47||ireal==50              %realization didn't worked out irate 4
%  	  if ireal==5 || ireal==21 || ireal== 29 || ireal==30 || ireal==32 || ireal==50
%     if ireal==21
        continue
    else
    
        %get CBTC
        filename = [pathCBTCs,basename,'_real',int2str(ireal),'_d',int2str(idpth),'_r',int2str(irate),'.dat'];
%         filename = [pathCBTCs,basename,'_real',int2str(ireal),'_d',int2str(idpth),'_r',int2str(irate),'_w',int2str(iwllp),'.dat'];
        fid      = fopen(filename);
        fgets(fid); fgets(fid);

        iwll  = 0;
        while 1 %~feof(fid)
            data = fgets(fid);
            if ischar(data)
                testZone = sscanf(data, '%c', 4);
                if strcmp(testZone,'ZONE') == 1 %change TP zone => next well
                    iwll = iwll+1;
                    k = 0;
                else
                    dataspl = sscanf(data, '%f%f%d%*f%*f%*f', 131); %either data or blank line
                    if size(dataspl) == [3 1]  %if data is [t,c,id]
                        k = k+1;
                        CBTC(k,1,iwll) = dataspl(1);
                        CBTC(k,2,iwll) = dataspl(2);
                    end
                end
            else
                break
            end
        end

        CBTC(:,1,:)=CBTC(:,1,:)./365;           %t in years
        CBTC(:,2,:)=CBTC(:,2,:)./(993270/100);  %mass in pourcentage of injected mass

        for iwll=1:nwell
            %plot CBTC
            ndata = nnz(CBTC(:,1,iwll));
            plott = CBTC(1:ndata,1,iwll);
            plotc = CBTC(1:ndata,2,iwll);

            %figure((nwell+2)*(irate-1)+1)
%             figure(20)
%             if ireal==1
%                 clf
%             end
%             semilogy(plott,plotc,'Color',[0.5 0.5 0.5])
%             hold on

            %analyse arrival times
            it=it+1;
            maxC = max(plotc);
%             for j=1:ndata
%                 for itana=1:size(tana,2)
%                     if plotc(j,1) < (tana(itana)/100)*maxC
%                         ti(it,itana) = plott(j,1);
%                     end
%                 end
%             end
            for itana=1:size(tana,2)
                for j=1:ndata
                    if plotc(j,1) > (tana(itana)/100)*maxC
                        if j>1
                            ti(it,itana) = plott(j-1,1)+((tana(itana)/100)*maxC-plotc(j-1,1))*(plott(j,1)-plott(j-1,1))/(plotc(j,1)-plotc(j-1,1));
                            break
                        else
                            ti(it,itana) = ((tana(itana)/100)*maxC)*(plott(j,1))/(plotc(j,1));
                            break
                        end
                    end
                end
            end
            
            for iw=1:size(win_min_ind,2)
                aw(it,iw) = ti(it,win_max_ind(iw))-ti(it,win_min_ind(iw));
            end
            
        end
        clear CBTC

    end
end
toc


%--------------------------------------------------------------------------
%PLOTS

pdfs_x = zeros(101,size(tana,2));
pdfs_y = zeros(101,size(tana,2));
% av_ti = zeros(size(tana,2),1);
% cv_ti = zeros(size(tana,2),1);
% var_ti = zeros(size(tana,2),1);

%get pdfs of ti
% for itana=1:size(tana,2)
%     
%     av_ti(itana,1) = mean(ti(:,itana));
%     cv_ti(itana,1) = std(ti(:,itana))/av_ti(itana,1);
%     var_ti(itana,1) = std(ti(:,itana))*std(ti(:,itana));
%     
%     pd = fitdist(ti(:,itana),'Kernel');
%     
%     nval = nnz(ti(:,itana));
%     mint = min(ti(1:nval,itana));
%     maxt = max(ti(1:nval,itana));
%     
%     x = mint:(maxt-mint)/100:maxt;
%     y = pdf(pd,x);
%     
%     %save pdf
%     pdfs_x(:,itana) = x;
%     pdfs_y(:,itana) = y;
% %     
%     %figure((nwell+2)*(irate-1)+2)
%     figure(21)
%     if itana==1
%         clf
%     end
%     plot(x,y,'Color',[1-itana/size(tana,2) 1-itana/size(tana,2) 1-itana/size(tana,2)],'LineWidth',2)
%     %plot(x,y,'Color','k','LineWidth',2)
%     %axis([0 400 0 0.014])
%     set(gca,'FontSize',18)
%     xlabel('arrival time [y]','FontSize',18);
%     ylabel('pdf','FontSize',18);
%     
%     hold on
% end

% get pdfs of aw
for iw=1:size(win_min_ind,2)
    
    pd = fitdist(aw(:,iw),'Kernel');
    
    nval = nnz(aw(:,iw));
    mint = min(aw(1:nval,iw));
    maxt = max(aw(1:nval,iw));
    
    x = mint:(maxt-mint)/100:maxt;
    y = pdf(pd,x);
    
    %save pdf
    pdfs_x(:,iw) = x;
    pdfs_y(:,iw) = y;
  
    %EXPORT TO FILE
    pathfile = ['D:\Box\batch_new\Plots\Outputs\pdfs\arriv_win\win_',int2str(tana(win_min_ind(iw))),'-',int2str(tana(win_max_ind(iw))),'\'];
    filename = [pathfile,'pdfArrWin','_d',int2str(idpth),'_r',int2str(irate),'.xlsx'];
    xlswrite(filename,pdfs_x(:,iw),'time')
    xlswrite(filename,pdfs_y(:,iw),'pdf')
end


%--------------------------------------------------------------------------
% pathfile = 'D:\Box\batch_new\Plots\Outputs\pdfs\';
% filename = [pathfile,'statTime_bw_qwei','_d',int2str(idpth),'_r',int2str(irate),'_w',int2str(iwllp),'.xlsx'];
% xlswrite(filename,pdfs_x,'time')
% xlswrite(filename,pdfs_y,'pdf')

fclose('all');

% %t vs wll_th
% wllinfo = dlmread([pathwllinfo,wllinfonam,'_d',int2str(idpth),'_r',int2str(irate),'.txt']);
% textscan(fidcbtc,'%f%f%f%f%f%f','Headerlines',3,'CommentStyle','@');
% for itana=1:size(tana,2)
%     figure(15+itana)
%     scatter( ti(:,itana),wllinfo(1:147,7) ) 
% end
