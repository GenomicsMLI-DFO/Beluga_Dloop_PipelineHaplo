# Beluga (_Delphinapterus leucas_) pipeline for D-loop haploype assignment

**Current version : 1.0.0**  
Check [this page](https://github.com/GenomicsMLI-DFO/Beluga_Dloop_PipelineHaplo/releases) for current and previous versions of the pipeline

__Main authors:__   Luca Montana, Geneviève Parent, and Benjamin Hornoy
__Affiliation:__    Fisheries and Oceans Canada (DFO)  
__Group:__          Laboratoire de génomique (Direction des sciences démersales et benthiques)  
__Location:__       Institut Maurice-Lamontagne  
__Contact:__        email: genevieve.parent@dfo-mpo.gc.ca  



- [Objective](#objective)
- [Summary](#summary)
- [Status](#status)
- [Contents](#contents)
  + [00_Data](#00_Data)
  + [01_Codes](#01_Codes)
  + [02_Results](#02_Results)
- [Methods](#methods)
- [Requirements](#requirements)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
Pipeline for generating haplotype libraries for short (234 bp) and long (615 bp) sequences of the mtDNA control region (D-loop), and provide haplotype assignment for beluga (*Delphinapterus leucas*) specimens


## Summary
Beluga from the Hudson Bay-Hudson Strait Complex are hunted by Nunavik and Nunavut communities for subsistence [1]. Within the Hudson Bay-Hudson Strait Complex, genetic studies highlighted distinct genetic groups, which are segregated during summer, using the mtDNA D-loop. Early studies using short mtDNA D-loop haplotypes (ca. 234 nucleotides) identified that eastern Hudson Bay and western Hudson Bay animals had distinct haplotype compositions [2-5]. More recently, longer haplotypes (615 nucleotides) highlighted the genetic distinctiveness of summering individuals from James Bay, Cumberland Sound, and Balcher Islands [6-10]. Because a mix of beluga from different summering areas are hunted by Inuit communities during spring and fall migration, as well as in winter, in the Hudson Strait and Ungava bay, it is imperative to correctly discriminate the summering area of origin of harvested belugas to understand the the proportion of animals summering in the eastern Hudson Bay [10-11].  
  
The R scripts included in this project guide users to generate consistent and reproducible haplotype libraries for the short (234 nucleotides) and long (615 nucleotides) haplotypes from available beluga mtDNA D-loop sequences, and consequently attribute the appropriate haplotype to each specimen based on their mtDNA D-loop sequence. This R pipeline contains data preparation, sequences alignment, detection of minimal sequence, haplotype library compilation, and haplotype assignation.


## Status
Under review


## Content
There are three folders (00_Data, 01_Codes, 02_Results) in this GitHub repository.

### 00_Data
**Home directory**  
* This pipeline is set to run on *genyoda*. Home directory is where the Rproj file is placed (*Beluga_Dloop_PipelineHaplo.Rproj*). From personal directory in *genyoda*, the path to access the project is: *../../media/genyoda/Extra_Storage/Projects/MOBELS/Beluga_Dloop_PipelineHaplo*
  
**/01_fasta**
* Complete mtDNA control region sequence (NCBI ID: *U18117.1*)  
* 615 nt consensus sequence created with Geneious  
  
**/02_dloop_clean**  
* Raw fasta files are provided here for n beluga.  
  
### 01_Codes
**0_DataPreparation.R**  
* Prepares sequences for haplotype libraries (short and long)  
* Mulitple sequences alignment  
* Cut 'complete' sequences to produce short (234 nt) and long (615 nt)  
    - Short: use specific position of the complete mtDNA D-loop region [12]  
    - Long: use consensus sequence to find cutting position  
* Fasta files in *project_name*/fasta
* Outputs will be placed in the home directory (*project_name*)  
  
**1_SNPsMinimalSequence.R**  
* Estimates number of polymorphisms, SNPs positions, limits of minimal sequence  
* Use minimal sequence to define haplotype in next scripts  
* Output will be placed in the home directory  

**2a_HaploLibrary_234.R** & **2b_HaploLibrary_615.R**
* Defines usable sequences
* Before executing, verify limits of minimal sequence to define haplotypes
* Output will be placed in home directory (*project_name*)  
* Compiles short haplotypes (234 nt) library from scratch (first time), or by extending an already existing one
* Verifies if sequences in library has ambiguities or empty sites - out of minimal sequence
* Creates 'libraries' directory that will be used to place new/updated haplotype libraries
* Libraries will be placed in *project_name*/libraries  
  
**3_AssignHaplo.R**  
* Assigns (short and long) haplotypes to each specimen 
* Uploads 'clean' sequences and libraries  
* Uses info on minimal sequence to assign haplotype to each sequence (for loop)  
* Creates new 'D-Loop' datasheet to be uploaded (manually) in ACCESS dataset
* Output will be placed in home directory  

### 02_Results
**00_Libraries**
* Lists reference librairies information for short and long haplotypes

**01_poly_seq_min**
* Provides the positions of SNPs considered in each reference library
 
**02_ACCESS**
* Provides the complete access information for all specimens at the Institut Maurice-Lamontagne (some without haplotypes)
  

## Methods
Please see the reference associated to this ReadMe below.


## Requirements
Download latest version of the following ACCESS data sheets: 'Specimens', 'D-loop'.
For externals to DFO, the reference library will be updated annually. Downloading the latest library prior haplotype identifiers assignment would be a good practice.


## Acknowledgements
Please see the reference below.


## References
This GitHub ReadMe is associated to this publication under review: Parent, G.J.*, Montana, L.*, Bonnet, C.*, Parent, É., Sauvé, C., St-Pierre, A.P., Watt, C., Hammill, M. 2024. Genetic monitoring program for beluga (Delphinapterus leucas) harvested in the Nunavik and Nunavut (Belcher Islands) regions. Can. Tech. Rep. Fish. Aquat. Sci. 0000 : v + x p.

These publications may also be useful to complete the understanding of this GitHub ReadMe.
1. Parent, G.J., Mosnier, A., Montana, L., Cortial, G., St-Pierre, A.P., Bordeleau, X., Lesage, V., Watt, C., Postma, L., Hammill, M.O. 2022. Reexamining populations of beluga in the Hudson Bay-Strait Complex and assessing the impact on harvests in Nunavik and Sanikiluaq management units.  DFO Can. Sci. Advis. Sec. Res. Doc. 2022/XXX.
2. Mosnier, A., Hammill, M.O., Turgeon, S., and Postma, L. 2017. Updated analysis of genetic mixing among beluga stocks in the Nunavik marine region and Belcher Islands area: information for population models and harvest allocation. DFO Can. Sci. Advis. Sec. Res. Doc. 2017/016. v + 15 p
3. Lillie, W. R., Gladden, J. G. B., and Tretiak, D. N. 1996. Amplification and sequencing of control region mitochondrial DNA from the beluga whale, Delphinapterus leucas. Can. Techn. Rep. Fish. Aquat. Sci., 2080, 8pp.  

