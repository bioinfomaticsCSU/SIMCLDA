function [feature_vecs]=pca_energy(simMat,para)
%pca_energy extracting primary feature vectors based on energy strategy 
%   Usage:  [feature_vecs]=pca_energy(simMat,para)
%	Inputs:
%			simMat: kernel of lncRNAs or diseases
%			para:	percent of energy for extracting primary feature vectors
%
%	Outputs:
%			feature_vecs: primary feature vectors of lncRNAs or diseases

    pca_rank=0;
    singular_mat=svd(simMat);
    for i=1:rank(simMat)
        if sum(singular_mat(1:i))>para*sum(svd(simMat))
            pca_rank=i;
            break;
        end
    end   
    [U_,S,V]=svds(simMat,pca_rank); 
    feature_vecs=V;     
end