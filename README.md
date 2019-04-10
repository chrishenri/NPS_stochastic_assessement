# Stochastic Assessment for Non-Point Source Contamination of Heterogeneous Aquifer: Instructions for Inputs and Outputs
Christopher V. Henri
chenri@ucdavis.edu

This document provides useful instructions to generate the input files of the different software used to model transport from a nonpoint source into a heterogeneous aquifer within a stochastic framework. 
We refer to the manuscripts attached to this document if further information about the mathematical and conceptual background of the study is needed. 


## 1.	T-PROGS
The same T-PROGS model is used for all studies made in the project. In our study, we do not use any data to condition the geostatistical model, but T-PROGS is proposing this option. Two set of parameters have to be generated in order to run the program: MCMOD inputs and TSIM inputs. MCMOD and TSIM inputs were generated using the GUI available for the TPROGS software. T-PROGS manual provides the scientific background and the detailed instructions on how to generate the needed inputs. A set of example parameter files and the tsim executable in provided in the folder 1_TPROGS. The software tsim will generate a single file with the facies spatial distribution for all realizations (output file with extension “.asc”).  
The only post-process required here to proceed with our stochastic analysis is to split this output file into a file per realization and convert the facies index to a corresponding value of hydraulic conductivity. A Matlab script does this: 

<ul>
<li>tsim_to_Kmat.mat: 	Generates hydraulic conductivity fields to be later used to generate the input of the Modflow model
</li>
</ul>


## 2.	MODFLOW
Some inputs are common to all realizations. It should be created once only, depending on the model main characteristics, and be used for any simulation. These files will be copied in the sub-folder used to run Modflow-2000 (MF2K) during the Monte Carlo simulation. These files, common to all realizations are: BAS, DIS, GMG, MNWI, OC. We refer to the Modflow online manual for details about these packages. An example of each file is provided in the folder 2_MF2K.
Some other input files depend on the hydraulic conductivity field, such as the recharge rate, the well package and the K-weighted fluxes leaving the domain and are, therefore, realization dependent. Some parameters used as input to compute the flow field is also dependent on the land use. 
These provided Matlab scripts help generate the series (for all realization) of input files: 

<ul>
<li>get_landuse_cons.m: Generates a series of RCH files, which are MF2K input files with spatial variability of the recharge rate; and generates the spatial distribution of the particle density used for the transport modeling. 
In our study, the local recharge rate depends on the soil and the crop types. The soil type is defined as the first (top) layer of the hydraulic conductivity field previously generated. 
	The Matlab script generates a random spatial distribution of a series of crop over the domain, which is also used to generates the spatial distribution of the particle density used for the transport modeling. 
	The crop and soil type dependent recharge rate and contaminant mass flux is obtained by means of a series of Hydrus-1D simulations. We refer to the attached manuscript for more details. </li>
	
<li>KweightedFHB.m: Generates a series of FHB packages, which specify the spatial distribution of prescribed fluxes to be applied at the bottom of the domain in order to simulate non-represented extraction. The local flux is proportional to the local hydraulic conductivity. </li>

<li>mnw2_pack.m: Generates the MNW2 packages for each realization with different pumping rate, screen length, and top depth. 3 extraction wells are implemented in each simulation. The well location is selected in order to always have 10 ft of gravel/sand for each 100 gpm of pumping rate. If this is not doable at the given well location, the algorithm changes this location until the criteria is fulfilled. </li>

<li>get_nam.m: Generates the name file (specifies MF2K input and output files names) for each realization. </li>
</ul>

## 3.	RW3D
The transport is solved using the random-walk particle-tracking method. We use the code RW3D, which is provided (compiled executable and source code) in its last version. The code presents a high versatility in the processes to be simulated: linear sorption (retardation), linear reaction network (first-order decay), multi-rate mass transfer, non-linear bi-molecular reaction network. We refer to the attached readme file for details about the specification of required inputs. In our Monte Carlo framework, all inputs are realization dependent and are generated using the following Matlab script: 

<ul>
<li>get_nameRW3D.m: Generates the name file (specifies RW3D input and output files names) for each realization. </li>

<li>get_paramRW3D.m: Generates the main input file specifying all parameters for each realization. Among other, this file is calling the flow field previously generated by MF2K and the file specifying the spatial distribution of the original particle density generated by the script get_landuse_cons.m. </li>
</ul>

## 4.	OUTPUTS
The objective of our studies is to better understand the spatio-temporal behavior of a series of management metrics adapted to nonpoint source contamination under highly uncertain conditions. To do so, we designed a series of visual tools to easily assess stochastically assessed travel times, capture zone (or contributing area) spatial and temporal extension, and contaminant levels at extraction wells. 


## 5.	Adaptation to pesticide contamination

