%K-weighted flux boundary
clear all

nreal = 50;

homo = 0; %0 if K heterogeneous; 1 otherwise
KVhomo = 0.13; %vertical K if homogeneous

basenameK = 'KmatMF';
pathK = 'D:\Box\batch_new\0_MF2K_inputs\K-mat\';

basenameout = 'box';
% pathout = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\fhb_pack_new\';
pathout = 'D:\Transient_Recharge\SteadyState\HOMOGENEOUS_K\get_Keq\0_MF2K_inputs\fhb_pack_nopump\';

% path_rch = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\rch_SS_pack\';
path_rch = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\rch_SS_pack\';
basename_rch = 'box_rch';

ncol = 120;
nrow = 60;
nlay = 625;

dx = 160.0;
dy = 100.0;
nwll = 3; %number of wells
Qwll = 0.0; %flux out a single well
%rchrate = 9.5E-4; %recharge rate


Header1 = '  1  7080  0  1  0  0  0';
Header2 = '  41  1.0  0';
Header3 = '  0.000000';
Header4 = '  41  1.000000  0';

for ireal=1:nreal
    
    %get average recharge rate for realization ireal
    fileRCH = [path_rch,basename_rch,'_',int2str(ireal),'.rch'];
    %rch = dlmread(fileRCH,' ',[5 1 5 20]); 
    rch = dlmread(fileRCH,' ',5,0); 
    rchrate = mean(rch(rch~=0));
    Qout = rchrate*dx*dy ;
    
    
    %get last layer of K-field
    Kmatl = zeros(nrow,ncol);
    if homo==0
        fileKmat = [basenameK,'_',int2str(ireal),'.dat'];
        Kmat     = dlmread(fullfile(pathK,fileKmat) ); 
        Kmatl    = Kmat((nlay-1)*nrow+1:nlay*nrow,:);
    end
    if homo==1
        Kmatl(:) = KVhomo;
    end
    
    %get K-weighted flux
    sumK     = sum(sum(Kmatl));
    Qtot     = Qout*ncol*nrow - Qwll*nwll;
    Q        = -Kmatl/sumK*Qtot;
    
    %open file rch package
    fileFHB = [pathout,basenameout,'_',int2str(ireal),'.fhb'];
    fileID  = fopen(fileFHB,'w');
    fprintf(fileID,'%c',Header1); fprintf(fileID,'\n');
    fprintf(fileID,'%c',Header2); fprintf(fileID,'\n');
    fprintf(fileID,'%c',Header3); fprintf(fileID,'\n');
    fprintf(fileID,'%c',Header4); fprintf(fileID,'\n');
    
    for irow=1:nrow
        for icol=2:ncol-1
            fprintf(fileID,'%c','  ');
            fprintf(fileID,'%i',nlay);
            fprintf(fileID,'%c','  ');
            fprintf(fileID,'%i',irow);
            fprintf(fileID,'%c','  ');
            fprintf(fileID,'%i',icol);
            fprintf(fileID,'%c','  0');
            fprintf(fileID,'\n');
            
            %Q_str = sprintf('%14.6e',-Q(irow,icol));
            %Q_str = strrep(Q_str, 'e+0','e+00');
            %Q_str = strrep(Q_str, 'e-0','e-00');
            %fprintf(fileID,'%s',Q_str);
            
            fprintf(fileID,'%14.6e',Q(irow,icol));
            fprintf(fileID,'\n');
        end
    end
    
end
fclose('all');