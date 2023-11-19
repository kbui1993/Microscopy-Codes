This Matlab Toolbox provides a set of functions to perform segmentation on
scanning tunnel microsocpy images. The methods and algorithms are described
in the paper

-K. Bui, et al. "Segmentation of Scanning Tunnel Microscopy Images Using Variational
Methods and Empirical Wavelets", 2019

https://link.springer.com/epdf/10.1007/s10044-019-00824-0?author_access_token=X7Qeci9q7wd9AR1fM7tYVve4RwlQNchNByi7wbcMAY69s1zBW2Cf_hewHCFXRC4WmygDLaBwjCkbEXfx2AytC3MLWJNCQE_xG7y7D_F6iOz67sKtBgJvP2QNQf7VWGEXEOn5iADNmrdD-WIO4UaYGg%3D%3D

To cite this paper,
```
@Article{Bui2019,
author="Bui, Kevin
and Fauman, Jacob
and Kes, David
and Torres Mandiola, Leticia
and Ciomaga, Adina
and Salazar, Ricardo
and Bertozzi, Andrea L.
and Gilles, J{\'e}r{\^o}me
and Goronzy, Dominic P.
and Guttentag, Andrew I.
and Weiss, Paul S.",
title="Segmentation of scanning tunneling microscopy images using variational methods and empirical wavelets",
journal="Pattern Analysis and Applications",
year="2019",
month="May",
day="05",
issn="1433-755X",
doi="10.1007/s10044-019-00824-0",
url="https://doi.org/10.1007/s10044-019-00824-0"
}
```
This toolbox is freely distributed and can be used without any charges for
research purposes, but please cite our papers when publishing your works.
The code should not be used for commercial purposes. 

For any questions, please contact kbui1993@gmail.com.

==========================Outside Toolboxes===============================
In the development of this toolbox, we require the following toolboxes developed by people 
outside the research group:

-Elad's Pseudo-Polar FFT toolbox: http://cs.techion.ac.il/~elad/software/
-Gilles' EWT toolbox: https://www.mathworks.com/matlabcentral/fileexchange/42141-empirical-wavelet-transforms

For the convenience of the user, the outside toolboxes are included in the 
Utilities subfolder in the Texture Segmentation folder. We do not claim that 
we contribute to any programming of the outside toolboxes.

===========================================================================

===========================Organization====================================

This toolbox is organized as follows:

Microscopy Codes
 |
 |-Cartoon Segmentation: functions for cartoon segmentation
 |-Cartoon+Texture: functions for cartoon+texture decomposition
 |-Images: contains the images used in the paper
 |-Segmentation Framework: contains a script to run segmentation framework on an image
 |-Test: contains scripts to run basic tests on the images in Images folder
 |-Texture Segmentation: functions for texture segmentation
     |-Utilities: contains outside toolboxes

===========================================================================
     
=============================Utilization===================================
The Test folders contain scripts to reproduce the segmentation results of 
the images in our paper. The entire_cartoon_script.m and the
entire_texture_script.m produce the cartoon and texture segmentation results,
respectively, for all sixteen images in the Image folder. There is a single
version if the user wants to run the segmentation algorithm onto one image.
The single version could be modified if the user wants to run it on his or
her own image. On the other hand, we have segmentation.m in Segmentation
Framework folder where the user can specify the image and the parameters
in the script and run the whole framework algorithm. 
