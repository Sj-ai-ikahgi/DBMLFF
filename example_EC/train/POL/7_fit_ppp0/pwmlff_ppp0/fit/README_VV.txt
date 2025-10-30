This directory contains some new subroutine for VV. 
Mostly, it has the select_mm_VV_forceMPI.f, which is a VV selection not only 
use the energy, but also use the force, to do the selection. Thhis is developed
since in the previous fitting, we found that, even if the energy is well described by VV, after 
added the force, the selected VV is not so good. Now, the force is also used to select the 
VV. Once again, the selected VV is output in: OUT.VV_index.1, OUT.VV_index.2 etc. 
One can use this file, to do further VV fitting. 

We also have a small program: feat_dist.f. This is run to get a set of new xp(2,id,j,itype), which 
is the Gaussian position and width to generate new feature: 
fnew=exp(-(f(j,itype)-xp(1,id,j,itype))^2/xp(2,id,j,itype)^2), from the input PCA (shifted and scaled) 
feature f(j,itype). The selection of xp is based on the distribution of f(j,itype). However, I found
the benefit is marginal (for the system I tested). Previoulsy, xp is fixed. Nevertheless, it does provide
something for fine tunning. Note, it will output OUT.xp, which will be needed (read) by the other 
fitting programs and calculating programs. 

The following is the procedure to generate the VV. 

(1) Run >python3 main.py, with parameters.py has: isCalcFeat=True, isFitLinModel=True. 
    This will generate fread_dfeat directory. cd to this directory. 

(2) create input file: select_VV.input (one example is in this directory). It has the following format:
----------------------
    10                 ! nloop
    2000               ! mm_tot
    20                 ! nimage_jump
    0                  ! include3
    20, 4, 2.0, 0.001  ! ndim1,ndim2,width,expd
----------------------
   nloop: the iteration loop to select VV. Larger nloop, better, but it will be more expensive. 
          nloop=10 is reasonable. 
   mm_tot: the total number of additional new feature, for each atom type. 
          mm_tot= 1000 to 2000 is recommended
   nimage_jump: The jump in the image to be used for selecting VV. Larger jump, it will be faster, 
                but it might not be so good. Note, this code can be rather slow. If nimage_jump=1, or 2 
                are chosen, it can take a few hours or one day. 
           nimage_jump=2 to 10 might be reasonable. 
  include3: whether to include the f1*f2*f3 as the new feature. Note, this can be extremely slow, and
            it might not be so useful. 
            include3=0 might be a good choice. 
  ndim1,ndim2,width,expd: The number of new feature for 1D and 2D. 1D means: fnew (exponential, see above). 
                          2D means: fnew(1)*fnew(2). 
             ndim1,ndim2=20,4 might be reasonable. Note, one should not make ndim2 too big! It will be 
             rather slow.  width is for exp(2,id,j,itype). Larger width, wider is the Guassian. width=1 or 2 
             is good. expd: is a power to select position of Guassian: xp(1,id,j,itype) according to the 
             feature distribution P(f): P(f)**expd. So, if expd is close to 0, then a more uniform grid will 
             be used. If expd is large, an more variance grid will be used. We found that expd=0.001 is 
             reasonable! 
-----------------------------------------
  run: >feat_dist.r
  It will take the ndim1,ndim2,width,expd to generate xp for 1D and 2D feature
  It will generate OUT.xp, to be used later. 
(3) run: >mpirun -n 8 select_mm_VV_forceMPI.r
   It will generate OUT.VV_index.itype (itype=1,2,..), these are the information for the VV new features. 
   It will also have an loop.inter (the loop iteration), it tells you over the nloop, when the new VV feature
   increases, how does the lost function: Lost=dE^2*weight_E + dF^2*weight_F reduce with the iteration. 
   The first line is the original linear fitting. However, this does not contain all the images. It only 
   contains the nimage_tot/nimage_jump images. 
   It also produces: energyVV.pred.type, forceVV.pred.type. One can plot them to see how good they are. 
   Note: energyVV.pred0.type, forceVV.pred0.type are the original linear fitting results. 
   Note: in the VV selection, we only used the energy: weight_E, and force: weight_F to do the selection. 
   The total energy weight_E0 is not used. 

(4) The above fitting also generates: linear_VV_fitB.ntype. So strictly speaking, one can already use this
    parameter file to do VV calculations, like MD etc. However, if nimage_jump is not 1, then not all the 
    data are used. Thus, we can redo a fitting with:
    >  fit_VV_forceMM.r
    followed by:
    > calc_VV_forceMM.r
    This will regenerate linear_VV_fitB.ntype using all the data point. 


    
