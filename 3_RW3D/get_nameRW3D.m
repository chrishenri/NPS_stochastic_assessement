%get all name file for RW3D
clear all

nrealK = 50;
ndepth = 3;
nrate  = 2;

path_nam = 'D:\Transient_Recharge\SteadyState\0_RW3D_inputs\nam\';
basenameNAM = 'rw3d';

pathinp  = 'D:\Transient_Recharge\SteadyState\0_RW3D_inputs\inp\';
pathbtc  = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\btc\';
pathcbtc = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\cbtc\';
pathmom  = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\spacemom\';
pathdbg  = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\dbg\';
pathex   = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\exit_part\';
%pathevol = 'D:\Box\batch\3_OUTPUTS\part\';

%INP  = ('rw3d.inp');
%BTC  = ('rw3d_BTCs.dat');
%CBTC = ('rw3d_CBTCs.dat');
EVOL = ('rw3d_part_evol.dat');
PATH = ('rw3d_part_path.dat');
%CMOM = ('rw3d_cart_spat_mom.dat');
SMOM = ('rw3d_xyz_plane_mom.dat');
PLAN = ('rw3d_xyz_plane_location.dat');
INDE = ('rw3d_dil_index.dat');
RMOM = ('rw3d_radial_spat_mom.dat');
TMOM = ('rw3d_temp_mom.dat');
DISP = ('rw3d_plane_dispersivity.dat');
TIME = ('rw3d_time.dat');
VELO = ('rw3d_velocities_TECPLOT.dat');
DBTC = ('rw3d_dBTCs.dat');
%DBG  = ('rw3d.dbg');
%EXIT = ('rw3d_exit_part.dat');


for irealK=1:nrealK
    for idepth=1:ndepth
        for irate=2:nrate
            % nam file name
            fname  = [path_nam,basenameNAM,'_',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.nam'];
            fileID = fopen(fname,'w');

            %inputs realization dependent
            INP  = [pathinp,'rw3d_',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.inp'];
            %outputs realization dependent
            BTC  = [pathbtc,'btc_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.dat'];
            CBTC = [pathcbtc,'cbtc_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.dat'];
            CMOM = [pathmom,'spacemom_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.dat'];
            DBG  = [pathdbg,'dbg_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.dbg'];
            EXIT = [pathex,'exit_part_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.dat'];
            %EVOL = [pathevol,'rw3d_part_evol',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.dat'];
            %print
            fprintf(fileID,'%c',INP);  fprintf(fileID,'\n');
            fprintf(fileID,'%c',BTC);  fprintf(fileID,'\n');
            fprintf(fileID,'%c',CBTC); fprintf(fileID,'\n');
            fprintf(fileID,'%c',EVOL); fprintf(fileID,'\n');
            fprintf(fileID,'%c',PATH); fprintf(fileID,'\n');
            fprintf(fileID,'%c',CMOM); fprintf(fileID,'\n');
            fprintf(fileID,'%c',SMOM); fprintf(fileID,'\n');
            fprintf(fileID,'%c',PLAN); fprintf(fileID,'\n');
            fprintf(fileID,'%c',INDE); fprintf(fileID,'\n');
            fprintf(fileID,'%c',RMOM); fprintf(fileID,'\n');
            fprintf(fileID,'%c',TMOM); fprintf(fileID,'\n');
            fprintf(fileID,'%c',DISP); fprintf(fileID,'\n');
            fprintf(fileID,'%c',TIME); fprintf(fileID,'\n');
            fprintf(fileID,'%c',VELO); fprintf(fileID,'\n');
            fprintf(fileID,'%c',DBTC); fprintf(fileID,'\n');
            fprintf(fileID,'%c',DBG);  fprintf(fileID,'\n');
            fprintf(fileID,'%c',EXIT); fprintf(fileID,'\n');
        end 
    end
end
fclose('all');