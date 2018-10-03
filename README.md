# 1dofResponseAnalysis
grasshopperでNewmarkのβ法により単質点系の応答解析を行うコンポーネントです.  
The  component which analyze of the single DOF mass point by Newmark β method in grasshopper.  

![コンポーネント画像](https://github.com/hiro-n-rgkr/1dofResponseAnalysis/blob/master/1dofResponseAnalysis/images/howtouse.PNG)

# LICENSE
MITライセンスです。詳細はLICENSE.mdを見てください。  
This software is released under the MIT License, see LICENSE.md  

# Operating Environment  
WindowsのRhino6向けに作っています。  
Make for Rhino6 of Windows

# Information Concerning the Author
hiron_rgkr  
twitter: https://twitter.com/hiron_rgkr

# How to Use  
## Commentary Video  
[![](http://img.youtube.com/vi/4HIVYRezRaI/sddefault.jpg)](https://www.youtube.com/watch?v=4HIVYRezRaI&t=3s)
   
## List of Component  
![List of Component](https://github.com/hiro-n-rgkr/1dofResponseAnalysis/blob/master/1dofResponseAnalysis/images/ListOfComponent.jpg)
  
## Component Appearance  
![Component Appearance](https://github.com/hiro-n-rgkr/1dofResponseAnalysis/blob/master/1dofResponseAnalysis/images/ComponentAppearance.PNG)  

## Function  
+ Calc T
  + 質量と剛性から固有周期、固有振動数、固有角振動数を計算  
  + A peculiar period and a natural frequency and a peculiar angular frequency are calculated by mass and rigidity
+ MTtoK  
  + 質量と固有周期から質量を計算  
  + Calculate rigidity from mass and a peculiar period  
+ KTtoM  
  + 剛性と固有周期から質量を計算  
  + Calculate mass from rigidity and a peculiar period  
+ MakeSinWave  
  + 入力波としてサイン波を作成  
  + Make a sin wave as an input wave  
+ 1dofRA 
  + 入力された諸元に対して、単質点系の応答解析をNewmarkβで計算し、応答結果とエネルギーを出力  
  + For input parameter, calculate 1dof response analysis in Newmarkβ method and output result  
+ ModelView  
  + 1質点モデルをRhino上に出力  
  + Output a model on Rhino  
+ ResultVeiw
  + Rhino上に出力した質点を応答結果に応じて動かす  
  + Move the model that output on Rhino depending on a result  

# Development Status  
+ ~2018/09/24
  + ver 0.1.00  
    基本的な機能を揃え、ver 0.1.0としてリリース  
    Prepare a basic function and release it as ver 0.1.0  
