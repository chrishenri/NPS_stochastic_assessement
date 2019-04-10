%generates parameters file for RW3D
clear all

pathout  = 'D:\Transient_Recharge\SteadyState\0_RW3D_inputs\inp\';
basename = 'rw3d';

pathwellinfo = 'D:\Transient_Recharge\SteadyState\0_MF2K_inputs\mnw2_pack\';
namewellinfo = 'mc_scenarios';

pathcbb = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\cbb\';
% fidwll = fopen(filenamewell,'rt');
% datawll=textscan(fidwll,'%f%f%f%f%f%f%f','Headerlines',1);
% fclose(fidwll);

pathmnw2 = 'D:\Transient_Recharge\SteadyState\3_OUTPUTS_WeakSink\Qmnw\';

pathpartloc = 'D:\Transient_Recharge\SteadyState\0_RW3D_inputs\npart_loc\';
namepartloc = 'npart_loc';

nreal  = 50;
ndepth = 3;
nrate  = 2;

nwll   = 3;

dx = 160.0;
dy = 100.0;
dz = 0.4;

nx = 120;
ny = 60;
nz = 625;

nt = 1;
dt = 146000.0;

%write inp file
for idepth=1:ndepth
    for irate=2:nrate
        
        %get well info
        filenamewell = [pathwellinfo,namewellinfo,'_d',int2str(idepth),'_r',int2str(irate),'.txt'];
        fidwll       = fopen(filenamewell,'rt');
        datawll      = textscan(fidwll,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fclose(fidwll);
        
        for ireal=1:nreal

            %generate rw3d input file name
            fname  = [pathout,basename,'_',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'.inp'];
            fileID = fopen(fname,'w');

            %get xwll,ywll,rwll,zbotwll,ztopwll
            xwll = zeros(nwll,1);
            ywll = zeros(nwll,1);
            rwll = zeros(nwll,1);
            ztopwll = zeros(nwll,1);
            zbotwll = zeros(nwll,1);
            for iwll=1:nwll
                xwll(iwll,1)    = datawll{1,5}(ireal*nwll-(nwll-iwll),1)*dx-dx/2;
                ywll(iwll,1)    = ny*dy - (datawll{1,6}(ireal*nwll-(nwll-iwll),1)*dy-dy/2);
                rwll(iwll,1)    = 1.0;
                ztopwll(iwll,1) = nz*dz - datawll{1,3}(ireal*nwll-(nwll-iwll),1);
                zbotwll(iwll,1) = ztopwll(iwll,1) - datawll{1,7}(ireal*nwll-(nwll-iwll),1);
            end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        %                                                   _/_/_/_/  _/      _/  _/_/_/_/  _/_/_/ 
        %                                                  _/    _/  _/      _/        _/  _/    _/ 
        %                                                 _/_/_/    _/  _/  _/    _/_/_/  _/    _/ 
        %                                                _/    _/  _/  _/  _/        _/  _/    _/ 
        %                                               _/    _/  _/_/_/_/_/  _/_/_/_/  _/_/_/ 
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % 			INPUTs for RW3D_Rx: REACTIVE TRANSPORT CODE BASED ON THE RANDOM WALK METHOD
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
            for i=1:9
                fprintf(fileID,'\n');
            end
            fprintf(fileID,'0');        fprintf(fileID,'\n'); 	%   debugging
            fprintf(fileID,'1   0');	fprintf(fileID,'\n'); 	%   Number of Aqueous and Mineral species
            fprintf(fileID,'A');        fprintf(fileID,'\n');  	%   Name of Aqueous species
            fprintf(fileID,'\n');                               %   Name of Mineral species 
            fprintf(fileID,'146000.0');      fprintf(fileID,'\n');   %   Simulation Time
            for i=1:3
                fprintf(fileID,'\n');
            end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . GEOMETRY .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf(fileID,'120    60    625');   fprintf(fileID,'\n');   %   nx,ny,nz
            fprintf(fileID,'not_used 			160.0 	1 	0');   fprintf(fileID,'\n');
            fprintf(fileID,'not_used 			100.0 	1 	0');   fprintf(fileID,'\n');
            fprintf(fileID,'not_used 			0.4 	1 	0');   fprintf(fileID,'\n');
            fprintf(fileID,'not_used 			1.0		1	0');   fprintf(fileID,'\n');
            fprintf(fileID,'0 0 0 0 0 0');   fprintf(fileID,'\n');
            for i=1:3
                fprintf(fileID,'\n');
            end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . TIME DISCRETIZATION .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf(fileID,'optimum_dt');   fprintf(fileID,'\n'); 	%   Method (constant_move, constant_time, one_time_x_cell, constant_Pe, constant_Pe_or_Cu, optimum_dt)
            fprintf(fileID,'0.1  0.2  0.3  0.1  0.1  0.1');   fprintf(fileID,'\n'); 	%   Dt, Cu, Pe, DaKINETIC, DaDECAY, DaMMT
                for i=1:3
                    fprintf(fileID,'\n');
                end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . ADVECTION . 
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------        
            fprintf(fileID,'T');   fprintf(fileID,'\n'); 	%   logical flag (F=no dispersion; T=yes)
            fprintf(fileID,'Eulerian');   fprintf(fileID,'\n'); 	%   method (Eulerian,Exponential)
            fprintf(fileID,'%c',pathcbb,'fluxes_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate)); 	%   qx(nx+1,ny,nz): file,const,ivar,flag   (0 => const, 1 => use file, 2 => use mf2k binary cell_by_cell flux file)
            fprintf(fileID,'.cbb   	1.0		0 	2 ');   fprintf(fileID,'\n');
            fprintf(fileID,'%c',pathcbb,'fluxes_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate)); 	%   qy(nx,ny+1,nz): file,const,ivar,flag   (0 => const, 1 => use file, 2 => use mf2k binary cell_by_cell flux file)
            fprintf(fileID,'.cbb   	1.0		0 	2 ');   fprintf(fileID,'\n');
            fprintf(fileID,'%c',pathcbb,'fluxes_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate)); 	%   qz(nx,ny,nz+1): file,const,ivar,flag   (0 => const, 1 => use file, 2 => use mf2k binary cell_by_cell flux file)
            fprintf(fileID,'.cbb   	1.0		0 	2 ');   fprintf(fileID,'\n');
            fprintf(fileID,'not_used     		0.3 	1 	0 ');   fprintf(fileID,'\n'); 	%   poro(nx,ny,nz): file,const,ivar,flag   (0 => const, 1 => use file)
            fprintf(fileID,'%d',nt);   fprintf(fileID,'\n'); 	%   NPER: 		Number of velocity stress periods
                for it=1:nt
                    fprintf(fileID,'%f',dt);
                    fprintf(fileID,'  1  1.0000  SS');   fprintf(fileID,'\n'); 	%   PERLEN:		length period, NSTP:number time steps, TSMULT: multiplier, SS/TR
                end
            fprintf(fileID,'T');   fprintf(fileID,'\n');
                for i=1:3
                    fprintf(fileID,'\n');
                end	

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . DISPERSION .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------        
            fprintf(fileID,'F');   fprintf(fileID,'\n'); 	%   logical flag (F=no dispersion; T=yes)
            fprintf(fileID,'not_used       		8.0     1 	0');   fprintf(fileID,'\n'); 	%   aL:  	file,const,ivar,flag   (0 => const, else => use file)
            fprintf(fileID,'not_used       		0.8		1 	0');   fprintf(fileID,'\n'); 	%   aTH: 	file,const,ivar,flag   (0 => const, else => use file)
            fprintf(fileID,'not_used       		0.01 	1 	0');   fprintf(fileID,'\n'); 	%   aTV:	file,const,ivar,flag   (0 => const, else => use file)
            fprintf(fileID,'not_used       		0.0 	1 	0');   fprintf(fileID,'\n'); 	%   Dm:  	file,const,ivar,flag   (0 => const, else => use file)
            fprintf(fileID,'not_used       		0.0 	1 	0');   fprintf(fileID,'\n'); 	%   DmTH:  	file,const,ivar,flag   (0 => const, else => use file)
            fprintf(fileID,'not_used       		0.0 	1 	0');   fprintf(fileID,'\n'); 	%   DmTV:  	file,const,ivar,flag   (0 => const, else => use file)
            fprintf(fileID,'1.0  1.0  1.0  1.0  1.0');   fprintf(fileID,'\n'); 	%   MULTa   (for each aqueous species)
            fprintf(fileID,'1.0  1.0  1.0  1.0  1.0');   fprintf(fileID,'\n'); 	%   MULTD   (for each aqueous species)
                for i=1:3
                    fprintf(fileID,'\n');
                end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . MULTIRATE MASS TRANSFER .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf(fileID,'F');   fprintf(fileID,'\n'); 	%   logical flag (F=no mass transfer; T=yes)
            fprintf(fileID,'spherical_diffusion');   fprintf(fileID,'\n'); 	%   type of mass transfer model (linear_sorption ; layered_diffusion ; spherical_diffusion ; multirate)
            fprintf(fileID,'10');   fprintf(fileID,'\n'); 	%   number of immobile zones or term for the Multirate series (for layered, cylindrical or spherical diffusion)
            fprintf(fileID,'not_used			0.037	1	0');   fprintf(fileID,'\n'); 	%   poro im_zone1:	file,const,ivar,flag (0 => const, else => use file) |porosity for each zone
            fprintf(fileID,'not_used			0.3		1	0');   fprintf(fileID,'\n'); 	%   alpha' zone1:	file,const,ivar,flag (0 => const, else => use file) |firt-order mass transfer rate for each zone
                for i=1:5
                    fprintf(fileID,'\n');
                end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . CHEMICAL REACTIONS .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % ... sorption
        % ----------------------------------------
            fprintf(fileID,'F');   fprintf(fileID,'\n'); 	%   logical flag (F=no mass transfer; T=yes)
            fprintf(fileID,'linear');   fprintf(fileID,'\n'); 	%   type of sorption (linear)
            fprintf(fileID,'not_used 			2.9		1	0');   fprintf(fileID,'\n'); 	%   Rm1: file,const,ivar,flag (0=> const,else => use file)
            fprintf(fileID,'not_used 			13.08	1	0');   fprintf(fileID,'\n'); 	%   Rim spe1-zone1
                for i=1:3
                    fprintf(fileID,'\n');
                end

        % ----------------------------------------
        % ... first order decay network
        % ----------------------------------------
            fprintf(fileID,'F');   fprintf(fileID,'\n'); 	%   logical flag (F=no mass transfer; T=yes)
            fprintf(fileID,'1');   fprintf(fileID,'\n'); 	%   number of species within the decay network
            fprintf(fileID,'B');   fprintf(fileID,'\n'); 	%   name of species (order controls the network if serial)
            fprintf(fileID,'serial');   fprintf(fileID,'\n'); 	%   reaction arquitecture (serial; network)
            fprintf(fileID,'not_used 			0.0042	1	0');   fprintf(fileID,'\n'); 	%   k1: file,const,ivar, flag (0=> const,else => use file)
            fprintf(fileID,'not_used 			0.00	1	0');   fprintf(fileID,'\n'); 	%   y11: file,const,ivar, flag (0=> const,else => use file)
            fprintf(fileID,'not_used 			0.0042	1	0');   fprintf(fileID,'\n'); 	%   kim spe1-zone1		|first-order decay rate in immobile zones
                for i=1:3
                    fprintf(fileID,'\n');
                end

        % ----------------------------------------
        % ... bimolecular reaction network
        % ----------------------------------------
            fprintf(fileID,'F');   fprintf(fileID,'\n'); 	%   logical flag (F=no mass transfer; T=yes)
            fprintf(fileID,'0   0');   fprintf(fileID,'\n'); 	%   number of chemical reactions
                for i=1:3
                    fprintf(fileID,'\n');
                end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . CONTROL SURFACE .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf(fileID,'%d',nwll);   fprintf(fileID,'\n'); 	%   nwell
            %fprintf(fileID,'W1  4600.0  3080.0  1.0  149.0  50.0  1  T');
            for iwll=1:nwll
                fprintf(fileID,'W');fprintf(fileID,'%d',iwll);fprintf(fileID,'  ');     % name,xwell,ywell,rwell,ztop,zbot,flag(1=removed,0=pass thru),SaveBTC  
                fprintf(fileID,'%f  %f  %f  %f  %f',xwll(iwll,1),ywll(iwll,1),rwll(iwll,1),zbotwll(iwll,1)+0.001,ztopwll(iwll,1));
                fprintf(fileID,'  1  T');
                fprintf(fileID,'\n');
            end
            fprintf(fileID,'mnw2_package');   fprintf(fileID,'\n');
            fprintf(fileID,'%c',pathmnw2,'Qmnw_real',int2str(ireal),'_d',int2str(idepth),'_r',int2str(irate),'.80');   fprintf(fileID,'\n');
            fprintf(fileID,'0');   fprintf(fileID,'\n'); 	%   nplane => planes: dist,type,flag(1=removed,pass thru)
                for i=1:3
                    fprintf(fileID,'\n');
                end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . INJECTIONS .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf(fileID,'1');   fprintf(fileID,'\n'); 	%   Number of Injections
            fprintf(fileID,'cells_file_particle_number    DIRAC');   fprintf(fileID,'\n'); 	%   Injection:  Type, Function (DIRAC/GENERAL)
            fprintf(fileID,'0  1');   fprintf(fileID,'\n'); 	%   properties: pmass, zone, species
            fprintf(fileID,'%c',pathpartloc,namepartloc,'_',int2str(ireal),'_ss.dat');   fprintf(fileID,'\n'); 	%   idwn,jdwn,kdwn,iup,jup,kup    
            %fprintf(fileID,'%c',pathpartloc,namepartloc,'_',int2str(ireal),'.dat');   fprintf(fileID,'\n'); 	%   idwn,jdwn,kdwn,iup,jup,kup  
            fprintf(fileID,'0');   fprintf(fileID,'\n'); 	%   Time of injection
                for i=1:3
                    fprintf(fileID,'\n');
                end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . RECIRCULATIONS .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf(fileID,'F');   fprintf(fileID,'\n'); 	%   logical flag (F=no mass transfer; T=yes)
            fprintf(fileID,'2');   fprintf(fileID,'\n'); 	%   Number of connection
            fprintf(fileID,'[ D1 AND D3 AND D5 AND D7 AND D9 AND D11 AND D13 AND D15 --> D2 AND D4 AND D6 AND D8 AND D10 AND D12 AND D14 ]');   fprintf(fileID,'\n'); 	%
            fprintf(fileID,'recircul_tf_odd.dat  	1.0');   fprintf(fileID,'\n'); 	%   Time function: file,const
            fprintf(fileID,'[ D2 AND D4 AND D6 AND D8 AND D10 AND D12 AND D14 --> D1 AND D3 AND D5 AND D7 AND D9 AND D11 AND D13 AND D15 ]');   fprintf(fileID,'\n'); 	%
            fprintf(fileID,'recircul_tf_even.dat  	1.0');   fprintf(fileID,'\n'); 	%   Time function: file,const
                for i=1:3
                    fprintf(fileID,'\n');
                end

        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        % . POST-PROCESSING AND OUTPUT OPTIONS .
        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf(fileID,'\n');            
            fprintf(fileID,'  1');   fprintf(fileID,'\n'); 	%   flag: Print Cartesian Spatial Moments at Snapshots	
            fprintf(fileID,'  0');   fprintf(fileID,'\n'); 	%   flag: Print Particle Cloud at Snapshots
            fprintf(fileID,'  146000.0   150   1.0');   fprintf(fileID,'\n'); 	%   times: OPTION_1: tlen,ntstep,tmult; OPTION_2: file name; OPTION_3: write "ALWAYS" to print every time step 
            fprintf(fileID,'\n');
            fprintf(fileID,'  0');   fprintf(fileID,'\n'); 	%   flag: Print Temporal Moments of BTCs     
            fprintf(fileID,'  1    100   plugin   -10.   0. 73000.');   fprintf(fileID,'\n'); 	%   flag: Print BTCs,ngrid,Kernel (BOX, TRIANGLE, GAUSS),bw (<0 => optimal),Min,Max
            fprintf(fileID,'  1    	1');   fprintf(fileID,'\n'); 	%   flag: Print CBTCs,frequency (particles/prints)
            fprintf(fileID,'  0    	1	-2');   fprintf(fileID,'\n'); 	%   flag: Print Path, frequency (moves/prints), particle (if < 0 => all particles)  
            for i=1:3
            	fprintf(fileID,'\n');
            end
        end
    end
end
fclose('all');