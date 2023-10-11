# Beluga (_Delphinapterus leucas_) pipeline for D-loop haploype assignment

**Current version : 1.0.0**  
Check [this page](https://github.com/GenomicsMLI-DFO/Beluga_Dloop_PipelineHaplo/releases) for current and previous versions of the pipeline

__Main authors:__   Benjamin Hornoy & Luca Montana  
__Affiliation:__    Fisheries and Oceans Canada (DFO)  
__Group:__          Laboratoire de génomique (Direction des sciences démersales et benthiques)  
__Location:__       Institut Maurice Lamontagne  
__Contact:__        email: luca.montana@dfo-mpo.gc.ca  



- [Objective](#objective)
- [Summary](#summary)
- [Status](#status)
- [Contents](#contents)
  + [Subsections within contents](#subsections-within-contents)
- [Methods](#methods)
  + [Subsections within methods](#subsections-within-methods)
- [Requirements](#requirements)
- [Caveats](#caveats)
- [Uncertainty](#uncertainty)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
Pipeline for generating haplotype libraries for short (234 bp) and long (615 bp) sequences of the mtDNA control region (D-loop), and provide haplotype assignment for beluga (*Delphinapterus leucas*) specimens


## Summary

Beluga from the Hudson Bay-Hudson Strait Complex are hunted by Nunavik and Nunavut communities for subsistence [1]. Within the Hudson Bay-Hudson Strait Complex, genetic studies highlighted distinct genetic groups, which are segregated during summer, using the mtDNA D-loop. Early studies using short mtDNA D-loop haplotypes (ca. 234 nucleotides) identified that eastern Hudson Bay and western Hudson Bay animals had distinct haplotype compositions [2-5]. More recently, longer haplotypes (615 nucleotides) highlighted the genetic distinctiveness of summering individuals from James Bay, Cumberland Sound, and Balcher Islands [6-10]. Because a mix of beluga from different summering areas are hunted by Inuit communities during spring and fall migration, as well as in winter, in the Hudson Strait and Ungava bay, it is imperative to correctly discriminate the summering area of origin of harvested belugas to understand the the proportion of animals summering in the eastern Hudson Bay [10-11].  
  
The R scripts included in this project guide users to generate consistent and reproducible haplotype libraries for the short (234 nucleotides) and long (615 nucleotides) haplotypes from available beluga mtDNA D-loop sequences, and consequently attribute the appropriate haplotype to each specimen based on their mtDNA D-loop sequence. This R pipeline contains data preparation, sequences alignment, detection of minimal sequence, haplotype library compilation, and haplotype assignation. To 


## Status
Ongoing-improvements


## Contents
Describe the contents of the repository. Are there multiple scripts or directories? What are there purpose and how do they relate to each other?

### Directories
  
**Home directory**  
  
* This pipeline is set to run on *genyoda*. Home directory is where the Rproj file is placed (*Beluga_Dloop_PipelineHaplo.Rproj*). From personal directory in *genyoda*, the path to access the project is: *../../media/genyoda/Extra_Storage/Projects/MOBELS/Beluga_Dloop_PipelineHaplo*
  

**/fasta**
  
* Complete mtDNA control region sequence (NCBI ID: *U18117.1*)  
* 615 nt consensus sequence created with Geneious  
  
**/libraries**  
  
* Short and long haplotype libraries are placed here  
* Both update libraries will be placed back here  
  

### Scripts

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
  


## Methods
What methods were used to achieve the purpose? This should be detailed as possible.  

### Subsections within methods
Often useful to organise the methods under multiple subheadings.


## Requirements
*Optional section.* List the input data requirements or software requirements to successfully execute the code.
Download latest version of the following ACCESS data sheets: 'Specimens', 'D-loop'

## Caveats
Anything other users should be aware of including gaps in the input or output data and warnings about appropriate use.


## Uncertainty
*Optional section.* Is there some uncertainty associated with the output? Assumptions that were made?


## Acknowledgements
*Optional section.* List any contributors and acknowledge relevant people or institutions


## References
1. Breton-Honeyman, K., Huntington, H.P., Basterfield, M., Campbell, K., Dicker, J., Gray, T., Jakobsen, A.E., Jean-Gagnon, F., Lee, D., Laing, R. and Loseto, L., 2021. Beluga whale stewardship and collaborative research practices among Indigenous peoples in the Arctic. Polar Research, 40.
2. Brown Gladden, J.G., Ferguson, M.M., and Clayton, J.W. 1997. Matriarchal genetic population structure of North American beluga whales _Delphinapterus leucas_ (Cetacea: Monodontidae). Mol. Ecol. 6: 1033-1046.
3. de March, B., Maiers, L. D, and Friesen, M. K. 2002. An overview of genetic relationships of Canadian and adjacent populations of belugas (_Delphinapterus leucas_) with emphasis on Baffin Bay and Canadian eastern Arctic populations. NAAMCO Sci. Publ. 4: 17-38
4. de March, B., Stern, G., and Innes, S. 2004. The combined use of organochlorine contaminant profiles and molecular genetics for stock discrimination of white whales (_Delphinapterus leucas_) hunted in three communities on southeast Baffin Island. J. Cetacean Res. Manage. 6: 241-250.
5. de March, B. G. E., and Postma, L. D. 2003. Molecular genetic stock discrimination of belugas (_Delphinapterus leucas_) hunted in eastern Hudson Bay, Northern Quebec, Hudson Strait, and Sanikiluaq (Belcher Islands), Canada, and comparisons to adjacent populations. Arctic 56: 111-124.
6. Turgeon, J. Duchesne, P. Postma, L.D., and Hammill, M.O. 2009. Spatiotemporal distribution of beluga stocks (_Delphinapterus leucas_) in and around Hudson Bay: Genetic mixture analysis based on mtDNA haplotypes. DFO Can. Sci. Advis. Sec. Res. Doc. 2009/011. iv + 14 p.
7. Turgeon, J., Duchesne, P., Colbeck, G.J.C., Postma, L., and Hammill, M.O. 2012. Spatiotemporal segregation among summer stocks of beluga (_Delphinapterus leucas_) despite nuclear gene flow: implication for an endangered population in eastern Hudson Bay (Canada). Conserv. Genet. 13(2): 419-433.
8. Postma, L.D., Petersen, S.D., Turgeon, J., Hammill, M.O., Lesage, V., and Doniol-Valcroze, T. 2012. Beluga whales in James Bay: a separate entity from eastern Hudson Bay belugas? DFO Can. Sci. Advis. Sec. Res. Doc. 2012/074. iii + 23 p.
9. Postma, L. D. 2017. Genetic diversity, population structure and phylogeography among belugas (_Delphinapterus leucas_) in Canadian waters: broad to fine-scale approaches to inform conservation and management strategies. PhD Thesis, University of Manitoba, 314 pp.
10. Parent, G.J., Mosnier, A., Montana, L., Cortial, G., St-Pierre, A.P., Bordeleau, X., Lesage, V., Watt, C., Postma, L., Hammill, M.O. 2022. Reexamining populations of beluga in the Hudson Bay-Strait Complex and assessing the impact on harvests in Nunavik and Sanikiluaq management units.  DFO Can. Sci. Advis. Sec. Res. Doc. 2022/XXX.
11. Mosnier, A., Hammill, M.O., Turgeon, S., and Postma, L. 2017. Updated analysis of genetic mixing among beluga stocks in the Nunavik marine region and Belcher Islands area: information for population models and harvest allocation. DFO Can. Sci. Advis. Sec. Res. Doc. 2017/016. v + 15 p
12. Lillie, W. R., Gladden, J. G. B., and Tretiak, D. N. 1996. Amplification and sequencing of control region mitochondrial DNA from the beluga whale, Delphinapterus leucas. Can. Techn. Rep. Fish. Aquat. Sci., 2080, 8pp.  

