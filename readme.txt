MBD.m is the algorithm for minimal branching decomposition of metabolic flux distributions described in:

Chan,S.H.J., Solem,C., Jensen,P.R. and Ji,P. (2014) Estimating biological elementary flux modes that decompose a flux distribution by the minimal branching property. Bioinformatics, doi: 10.1093/bioinformatics/btu529.

It aims to estimate more biological EFMs to decompose a flux distribution.
You need to supple a COBRA model (with 'S' and 'rev' at least) and a flux vector.
If an EFM matrix is not supplied, the function will try to identify EFMtool and calculate the EFM matrix.
(Terzer, M. and Stelling, J., "Large scale computation of elementary flux modes with bit pattern trees", Bioinformatics, July 28, 2008. doi:10.1093/bioinformatics/btn401) 
(freely available at http://www.csb.ethz.ch/tools/efmtool)

In 'example.mat', there are two structures:
 
1. 'Cardiomyocyte' is the structure containing the first 10
       optimal MBDs, 35 MinEFMDs and the problem structure that further 
       contains the rediced network, the flux distribution and the EFM matrix. 
       The original network data can be downloaded in the supplementary 
       information of the paper:
       Vo,T.D. and Palsson,B.O. (2006) Isotopomer analysis of myocardial 
       substrate metabolism: a systems biology approach. Biotechnol. Bioeng., 
       95, 972–83.
       (http://gcrg.ucsd.edu/InSilicoOrganisms/Cardiomyocyte)
       (http://onlinelibrary.wiley.com/doi/10.1002/bit.21063/suppinfo)
 
2. 'lactis' is the structure containing the first 10 optimal MBDs, 
       and the 10 MinEFMDs compared in the paper and the problem structure 
       that further contains the rediced network, the flux distribution 
       and the EFM matrix. The original network data can be downloaded 
       in the supplementary information of the paper:
       Flahaut,N. a L. et al. (2013) Genome-scale metabolic model for 
       Lactococcus lactis MG1363 and its application to the analysis of flavor
       formation. Appl. Microbiol. Biotechnol., 97, 8729–39.
       (http://link.springer.com/article/10.1007%2Fs00253-013-5140-2)


Please cite the paper if you has used the algorithm in your publication.
I earnestly hope it can help you. Please feel free to contact me if you have any question.
Thank you!

Siu Hung Joshua Chan
email: joshua.chan@connect.polyu.hk