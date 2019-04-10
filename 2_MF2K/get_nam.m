%get all name file for MF2K
clear all

nrealK = 50;
ndepth = 3;
nrate  = 2;

path_nam = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\nam_files\';
basenameNAM = 'mf2k_ss';

path_fhbpck = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\fhb_pack\';
basenameFHB = 'box';

path_rchpck = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\rch_SS_pack\';
basenameRCH = 'box_rch';

path_lst = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\lst\';
basenameLST = 'box';

path_cbb = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\cbb\';
basenameCBB = 'fluxes';


path_baspack = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\';
nameBAS = 'box.bas';

path_lpfpck  = 'D:\Box\batch_new\0_MF2K_inputs\lpf_pack\';
basenameLPF = 'box';

path_mnw2pck = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\mnw2_pack\';
basenameMNW = 'box';


GLOB  = ('GLOBAL 3 box.glo');
%LIST  = ('LIST 2 boxH2.lst');
%BAS   = ('BAS6 1 boxH2.bas');
BAS   = ['BAS6 1 ',path_baspack,nameBAS];
%LPF   = ('LPF 11 D:\Box\batch\1_MF2K\boxH2.lpf');
%FHB = ('FHB 41 box.fhb');
DIS   = ('DIS 29 box.dis');
%ZONE  = ('ZONE 40 boxH2.zone');
%WEL   = ('WEL 12 D:\Box\batch\1_MF2K\boxH2.wel');
%RCH   = ('RCH 18 box.rch');
OC    = ('OC 22 box.oc');
GMG   = ('GMG 19 box.gmg');
%MNW2 =  ('MNW2 51 box.mnw2');
MNWI = ('MNWI 61 box.mnwi');

%DATA1 = ('DATA(BINARY) 50 D:\Box\batch\1_MF2K\fluxes.cbb');
DATA2 = ('DATA(BINARY) 30 box.hds');
DATA3 = ('DATA(BINARY) 12 box.ddn');

for irealK=1:nrealK
    for idepth=1:ndepth
        for irate=2:nrate
            % nam file name
            fname  = [path_nam,basenameNAM,'_',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.nam'];
            fileID = fopen(fname,'w');

            %inputs realization dependent
            %GLOB = ['GLOBAL 3 ',path_glo,basenameGLO,'_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.glo'];
            LIST = ['LIST 2 ',path_lst,basenameLST,'_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.lst'];
            LPF = ['LPF 11 ',path_lpfpck,basenameLPF,'_real',int2str(irealK),'.lpf'];
            MNW = ['MNW2 51 ',path_mnw2pck,basenameMNW,'_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.mnw2'];
            FHB = ['FHB 41 ',path_fhbpck,basenameFHB,'_',int2str(irealK),'.fhb'];
            RCH = ['RCH 18 ',path_rchpck,basenameRCH,'_',int2str(irealK),'.rch'];
            %outputs realization dependent
            DATA1 = ['DATA(BINARY) 50 ',path_cbb,basenameCBB,'_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.cbb'];
            %DATA2 = ['DATA(BINARY) 30 ',path_hds,basenameHDS,'_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.hds'];
            %DATA3 = ['DATA(BINARY) 50 ',path_ddn,basenameDDN,'_real',int2str(irealK),'_d',int2str(idepth),'_r',int2str(irate),'.ddn'];
            %print
            fprintf(fileID,'%c',GLOB); fprintf(fileID,'\n');
            fprintf(fileID,'%c',LIST); fprintf(fileID,'\n');
            fprintf(fileID,'%c',BAS); fprintf(fileID,'\n');
            fprintf(fileID,'%c',LPF); fprintf(fileID,'\n');
            fprintf(fileID,'%c',DIS); fprintf(fileID,'\n');
            fprintf(fileID,'%c',FHB); fprintf(fileID,'\n');
            fprintf(fileID,'%c',RCH); fprintf(fileID,'\n');
            fprintf(fileID,'%c',OC); fprintf(fileID,'\n');
            fprintf(fileID,'%c',GMG); fprintf(fileID,'\n');
            fprintf(fileID,'%c',MNW); fprintf(fileID,'\n');
            fprintf(fileID,'%c',MNWI); fprintf(fileID,'\n');
            fprintf(fileID,'%c',DATA1); fprintf(fileID,'\n');
            fprintf(fileID,'%c',DATA2); fprintf(fileID,'\n');
            fprintf(fileID,'%c',DATA3); fprintf(fileID,'\n');
        end 
    end
end
fclose('all');