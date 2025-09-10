# Broughton physical classification for kelp version 2

__Created:__      2024/05/29 Cloned from DFO classification to support LSSM work.
__Main author:__  Edward Gregr  
__Affiliation:__  SciTech Environmental Consulting   
__Location:__     Vancouver, BC   
__Contact:__      e-mail: ed@scitechconsulting.com | tel: 604-612-8324
__Last Update:__  2025/09/03   
__Version:__      R version 4.3.2 x64

- [Objective](#objective)
- [Summary](#summary)

## Objective
Building on the k-means clustering of physical layers done for DFO/MSEA, this 
project focuses on building clusters for the Broughton region to support the 
creation of relevant oceanographic contexts. 

## Summary   
2024/05/29: Initially cloned from DFO work to allow clean delivery and further 
development. However, DFO work showed limitations of clustering without any 
obvious goals. Considerable time therefore spent on standardizing and 
normalizing the date, and RMD script development. 

2024/09/22: Returned to the Broughton clustering, transferring what was 
completed as part of the DFO classifications. Paused in Oct 2024 because of
distractions and unsatisfying results. 

2025/07/15: Return with intent to use FVCOM data directly, and to use learnings
from completed DFO work. Began with the direct loading of the Bianucci FVCOM 
data as I wanted some clear current data. First finding was that the SPECTRAL 
lab interpolations to 20 m were inappropriate, particularly across the the 
Village Sea, which contains no FVCOM points. 

2025/09/03: Had a good push in July with ChatGPT support to load FVCOM data from
NetCDF files, interpolation of currents from the edges to the nodes, and getting
some good classification results. Now need a plan to wrap this up. 


# fin.