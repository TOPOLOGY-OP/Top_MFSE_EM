# Top_MFSE_EM
The present code is supplementary to the corresponding paper:Anisotropic material-field series expansion for topological design  of optical metalens

1. The present code uses the MMA Matlab routines (which is widely used in the structural optimization community) as the optimizer. 
The present code calls the subroutine mmasub.m (“Version September 2007”) from the MMA Matlab routines, which in turn calls
another MMA subroutine subsolv.m (“Version Dec 2006”). 
As the authors of the MMA matlab routines require, the MMA Matlab routines can be obtained by contacting Prof. Svanberg (krille@math.kth.se). 

2. Once the MMA matlab routines are obtained, some modifications of the optimization parameters need to be done as follows to run the design problem in the present paper:
(i)  set the parameter asyinit = 0.1, 
(ii) change the move limits of the design variables by replacing the 
     original lines:

        if iter < 2.5 
          low = xval - asyinit*(xmax-xmin); 
          upp = xval + asyinit*(xmax-xmin); 
        else 
          zzz = (xval-xold1).*(xold1-xold2); 
          factor = eeen; 
          factor(find(zzz > 0)) = asyincr; 
          factor(find(zzz < 0)) = asydecr; 
          low = xval - factor.*(xold1 - low); 
          upp = xval + factor.*(upp - xold1); 
          lowmin = xval - 10*(xmax-xmin); 
          lowmax = xval - 0.01*(xmax-xmin); 
          uppmin = xval + 0.01*(xmax-xmin); 
          uppmax = xval + 10*(xmax-xmin); 
          low = max(low,lowmin); 
          low = min(low,lowmax); 
          upp = min(upp,uppmax); 
          upp = max(upp,uppmin); 
        end

     with

        if iter < 2.5 
          low = xval - asyinit*1; 
          upp = xval + asyinit*1; 
        else 
          zzz = (xval-xold1).*(xold1-xold2); 
          factor = eeen; 
          factor(find(zzz > 0)) = asyincr; 
          factor(find(zzz < 0)) = asydecr; 
          low = xval - min(2,factor.*(xold1 - low)); 
          upp = xval + min(2,factor.*(upp - xold1)); 
        end
