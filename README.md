# 1dofResponseAnalysis
The  component which analyze of the single DOF mass point by Newmark β method in grasshopper.  

![コンポーネント画像](https://github.com/hiro-n-rgkr/1dofResponseAnalysis/blob/master/1dofResponseAnalysis/images/howtouse.PNG)

# LICENSE
This software is released under the MIT License, see LICENSE.md  

# Operating Environment  
Make for Rhino6 of Windows

# Information Concerning the Author
hiron_rgkr  
twitter: https://twitter.com/hiron_rgkr

# How to Use     
## List of Component  
![List of Component](https://github.com/hiro-n-rgkr/1dofResponseAnalysis/blob/master/1dofResponseAnalysis/images/ListOfComponent.jpg)
  
## Component Appearance  
![Component Appearance](https://github.com/hiro-n-rgkr/1dofResponseAnalysis/blob/master/1dofResponseAnalysis/images/ComponentAppearance.PNG)  

## Function  
+ Calc T
  + A peculiar period and a natural frequency and a peculiar angular frequency are calculated by mass and rigidity
+ MTtoK  
  + Calculate rigidity from mass and a peculiar period  
+ KTtoM  
  + Calculate mass from rigidity and a peculiar period  
+ MakeSinWave  
  + Make a sin wave as an input wave  
+ 1dofRA 
  + For input parameter, calculate 1dof response analysis in Newmarkβ method and output result  
+ ModelView  
  + Output a model on Rhino  
+ ResultVeiw
  + Move the model that output on Rhino depending on a result  

# Development Status  
+ ~2018/09/24
  + ver 0.1.00  
    Prepare a basic function and release it as ver 0.1.0  
    
# Video  
[![](http://img.youtube.com/vi/4HIVYRezRaI/sddefault.jpg)](https://www.youtube.com/watch?v=4HIVYRezRaI&t=3s)
