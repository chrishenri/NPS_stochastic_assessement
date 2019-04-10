clear all
%from TSIM simulation (categories in gslib format) to K-matrix (in MODFLOW format)

%---------------------------------
%PARAMETERS, TO BE DEFINED
nreal = 1; %number of realization
pathascfile = 'D:\Box\batch_new\TPROGS\tsim\output\'; %path with TPROGS output (*.asc)
namefile = 'tsim_box.asc'; %name of the TPROGS output (*.asc)

%output file name
basenameSgems = 'lnK_boxSgems'; %name of the output file (in SGEMS format)
pathMF        = 'D:\Box\batch_new\0_MF2K_inputs\K-mat\'; %path of the output file (in SGEMS format)
basenameMF    = 'KmatMF'; %name of the output file (in Modflow format)

%hydrau. parameters
ncat = 4; %number of categories used in the TPROGS model

nx = 120; %number of cell in x
ny = 60;  %number of cell in y
nz = 625; %number of cell in z

K(1,1) = 200.0 ; %Hydraulic conductivity values for gravel
K(2,1) = 50.0 ; %Hydraulic conductivity values for sand
K(3,1) = 0.5 ; %Hydraulic conductivity values for muddy sand
K(4,1) = 0.01 ; %Hydraulic conductivity values for mud

PRINTK_GSLIB = 0 ; %1 if print output in GSLIB format; 0 otherwise
PRINTK_MF    = 1 ; %1 if print output in Modflow format; 0 otherwise

%END OF PARAMETERS
%---------------------------------


for ireal=1:nreal

    %input filename
    namereal = [pathascfile,namefile,int2str(ireal)];
    fileID   = fopen(namereal,'r'); %input file (from TSIM)

    %define output files name
    if PRINTK_GSLIB == 1
        nameSgems = [basenameSgems,'_',int2str(ireal),'.dat'];
        fileIDK = fopen(nameSgems,'a');
    end

    if PRINTK_MF == 1
        nameMF = [pathMF,basenameMF,'_',int2str(ireal),'.dat'];
        fileOut = nameMF;
    end
    %--------------------------------
    ndata = nx*ny*nz;

    data = fscanf(fileID,'%u');

    Kmat   = zeros(nx,ny,nz); 
    KmatMF = zeros(ny,nx,nz); 

    i=4;
    for iz=1:nz
        for iy=1:ny
            for ix=1:nx
                Kmat(ix,iy,iz) = K(data(i,1),1); % GSLIB format
                if PRINTK_GSLIB == 1
                    fprintf(fileIDK,'%4.2e\n', log(Kmat(ix,iy,iz)));
                end
                i=i+1;
            end
        end
    end
    if PRINTK_GSLIB == 1
        fclose(fileIDK);
    end
        
    if PRINTK_MF == 1
        %reorganize the K-matrix
        for ix=1:nx
            for iy=1:ny
                for iz=1:nz
                    KmatMF(iy,ix,iz) = Kmat(ix,ny-iy+1,nz-iz+1); % MF format
                end
            end
        end

        %write the K-matrix to be inputted in Modflow (y and z are upside-down)
        for iz=1:nz
            KmatPrint(:,:) = KmatMF(:,:,iz);
            dlmwrite(fileOut,KmatPrint,'delimiter','\t','-append')
            fid = fopen(fileOut, 'a');
            fprintf(fid, '\n');
            fclose(fid);
        end
    end

end
fclose('all');