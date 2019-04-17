% Run batch in parallel framework

% To do before running the batch: 
% .. for nreal K-field:
% 1. generate nreal asc files (tsim)
% 2. generate nreal K-fields (tsim_to_Kmat.m)
% 3. generate nreal lpf packages (lpf_pack.m)
% .... for nrate, for ndepth
% 4. generate nreal*nrate*ndepth well packages (well_pack.m)
% 5. generate nreal*nrate*ndepth name files (get_nam.m)
% 6. generate RW3D input parameter files (get_paramRW3D.m)
% 7. generate RW3D name files (get_nameRW3D)

% What the batch is doing?
% 1. create batch file for each realization that
%   a. copy MF2K nam file
%   b. runs MF2K
%   c. copy cbb file into RW3D folder
%   d. copy RW3D parameter file
%   e. runs RW3D
%   f. copy CBTC, BTC, spatial moments in corresponding folder
% 2. run batch


clear all
clc

pathbat      = ('D:\Transient_Recharge\SteadyState\run_batch\batches\');
pathMF2k     = ('D:\Transient_Recharge\SteadyState\1_MF2K\');
pathMF2k_inp = ('D:\Transient_Recharge\SteadyState\0_MF2K_inputs\');
pathRW       = ('D:\Transient_Recharge\SteadyState\2_RW3D\');

% pathDBG = ('D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\dbg\');

nreal  = 50;
ndepth = 3;
nrate  = 2;

mycluster = parcluster('local');
mycluster.NumWorkers = 12;

parpool('local',10);

tic
parfor ireal=1:nreal
% for ireal=1:nreal   
    %ireal
    
    for idepth=3:ndepth
        for irate=2:nrate
            
%             dbgfile = [pathDBG,'dbg_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'.dbg'];
%             if ~exist(dbgfile, 'file')

                %Create folders - copy files that are not realization dependent
                %.. MK2k
                foldMF = [pathMF2k,'MF2k_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'\'];
                if ~exist(foldMF, 'dir')
                    mkdir(foldMF);
                end

                if ~exist([foldMF,'MF2k.exe'], 'file')
                    copyfile([pathMF2k,'MF2k.exe'],foldMF)
                end

                if ~exist([foldMF,'box.dis'], 'file')
                    copyfile([pathMF2k_inp,'box.dis'],foldMF)
                end

%                 if ~exist([foldMF,'box.rch'], 'file')
%                    copyfile([pathMF2k,'box.rch'],foldMF)
%                 end

                if ~exist([foldMF,'box.oc'], 'file')
                    copyfile([pathMF2k_inp,'box.oc'],foldMF)
                end

                if ~exist([foldMF,'box.gmg'], 'file')
                    copyfile([pathMF2k_inp,'box.gmg'],foldMF)
                end

                if ~exist([foldMF,'box.mnwi'], 'file')
                    copyfile([pathMF2k_inp,'box.mnwi'],foldMF)
                end

                %.. RW3D
                foldRW = [pathRW,'RW3D_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'\'];
                if ~exist(foldRW, 'dir')
                    mkdir(foldRW);
                end

                if ~exist([foldRW,'RW3D_Rx.exe'], 'file')
                    copyfile([pathRW,'RW3D_Rx.exe'],foldRW)
                end

                %---------------------------------------------------
                %Create batch
                batname = [pathbat,'batch_',int2str(ireal),'.bat'];
                fileID  = fopen(batname,'w');
                fprintf(fileID,'%c','@echo off');                                                   fprintf(fileID,'\n');

                %..set pathes
                fprintf(fileID,'%c','if "%1" == "}{" goto %2');                                     fprintf(fileID,'\n');
                fprintf(fileID,'%c','set PathBatch="D:\Transient_Recharge\SteadyState\run_batch\batches"');                        fprintf(fileID,'\n');
                fprintf(fileID,'%c','set PathMF2K="D:\Transient_Recharge\SteadyState\1_MF2K\MF2k_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'"'); fprintf(fileID,'\n');
                fprintf(fileID,'%c','set PathNAM="D:\Transient_Recharge\SteadyState\0_MF2K_inputs\nam_files"');          fprintf(fileID,'\n');
                fprintf(fileID,'%c','set PathQmnw="D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\Qmnw"');          fprintf(fileID,'\n');
                fprintf(fileID,'%c','set PathRWNAM="D:\Transient_Recharge\SteadyState\0_RW3D_inputs\nam"');              fprintf(fileID,'\n');
                fprintf(fileID,'%c','set PathRW3D="D:\Transient_Recharge\SteadyState\2_RW3D\RW3D_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'"');	fprintf(fileID,'\n');

                %..batch
                fprintf(fileID,'%c',':----------- start programs to run in a batch here');          fprintf(fileID,'\n');
                fprintf(fileID,'\n');
                fprintf(fileID,'%c',':***** run MF2K');                                             fprintf(fileID,'\n');
                fprintf(fileID,'%c','cd %PathMF2K%');                                               fprintf(fileID,'\n');
                fprintf(fileID,'%c','MF2k %PathNAM%\mf2k_ss_',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'.nam');    fprintf(fileID,'\n');
                fprintf(fileID,'%c','xcopy fort.80 %PathQmnw%\Qmnw_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'.80* /Y');  fprintf(fileID,'\n');
                fprintf(fileID,'\n');
                fprintf(fileID,'\n');
                fprintf(fileID,'%c',':***** Run RW3D and get CBTCs');                               fprintf(fileID,'\n');
                fprintf(fileID,'%c','cd %PathRW3D%');                                               fprintf(fileID,'\n');
                fprintf(fileID,'%c','RW3D_Rx %PathRWNAM%\rw3d_',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'.nam');  fprintf(fileID,'\n');
                fprintf(fileID,'\n');
                fprintf(fileID,'\n');
                fprintf(fileID,'%c','cd %PathBatch%');                                              fprintf(fileID,'\n');
                fprintf(fileID,':-----------------------------------------------');             	fprintf(fileID,'\n');

                %Run batch
                % status = dos(batname);
                [status, temp] = dos(batname);

                rmdir(foldMF,'s');
                rmdir(foldRW,'s');
%             end
            
        end
    end
end
fclose('all');
toc