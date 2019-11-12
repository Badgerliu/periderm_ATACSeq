library(gkmSVM)


genNullSeqs_drerio('fish_training_no_promoters.bed', nMaxTrials = 10, xfold = 1, genomeVersion='danRer7', ouputPosFastaFN = 'fish_training_np_pos.fa', outputBedFN = 'fish_training_np_neg_1x.bed', outputNegFastaFN = 'fish_training_np_neg_1x.fa')
gkmsvm_kernel('fish_training_np_pos.fa', 'fish_training_np_neg_1x.fa', 'fish_training_np_1x_kernel.out')
gkmsvm_trainCV('fish_training_np_1x_kernel.out', 'fish_training_np_pos.fa', 'fish_training_np_neg_1x.fa', svmfnprfx = 'fish_training_np_1x', outputCVpredfn = 'fish_training_np_1x_cvpred.out', outputROCfn = 'fish_training_np_1x_roc.out')
gkmsvm_classify('nr10mers.fa', svmfnprfx = 'fish_training_np_1x', 'fish_training_np_1x_weights.out')
