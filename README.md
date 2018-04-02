License
=========
Copyright (C) 2017 Jianxin Wang(jxwang@mail.csu.edu.cn),Chengqian Lu(chengqlu@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn),Chengqian Lu(chengqlu@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083

Type: Package
Title: Prediction of lncRNA-disease associations based on inductive matrix completion
=================
Description: This package implements the SIMCLDA algorithm with inductive matrix completion framework, predicting lncRNA-disease 
associations.

Files:
1.Dataset

1) lncSim.mat and disSim_Jaccard.mat store lncRNA similarity matrix and disease similarity matrix, respectively;

2) interMatrix.mat stores known lncRNA-disease association information;

3) lncRNA_Name.txt and diseases_Name.txt store lncRNA ids and disease ids, respectively;

2.Code
1) gKernel.m: function computing Gaussian interaction profile kernel;

2) pca_energy.m: function extracting feature vectors via PCA;

3) SIMC.m : function completing matrix;

4) SIMCLDA: predict potential lncRNA-disease associations; 

3.