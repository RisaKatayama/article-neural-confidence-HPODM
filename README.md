# article-neural-mechanisms-prior-uncertainty-hierarchical-decision

This repository contains data and code supporting the following paper:

***

Unthresholded group-level statistical maps supporting our results are available on NeuroVault: [https://neurovault.org/collections/NZJMDMFQ/](https://neurovault.org/collections/FUAPHGRH/).

#### Folders
* _data_
  
  contains anonymised behavioural data for the behavioral experiment.
  * _behavioral_ contains data files from the behavioral experiments (outside of the scanner)
  * _scanning_ contains data files from the fMRI experiments
  * bhvdata_**.RData contains summarized data from the behavioral experiments (_bhv_) or the fMRI experiments (_mri_) respectively

* *res_model*

  contains fitting results using each type of proposed/alternative models

* *timeseries_trl*

  contains BOLD timeseries (0-1 scaled and mean-centered) in regions of interets for each participants during [-2s, 18s] from the decision onset (data matrix: _N_ trial x _T_ timepoints)

  s*X*_*ROI*_8mm.mat, for data in *ROI* (left ventral anterior insular cortex [vAIC] or dorsal anterior cingulate cortex [dACC]: please refer to *Supplementary Table 2*) for Subject _X_

* *PPI*

   contains psychophysiological interaction timeseries in regions of interets for each participants

   _s*X*_ folder includes data from Subject _X_

* *codes*

  contains scripts designed to reproduce model-fitting results (_Model_fitting.R_: please refer to **Methods**) and each figure (_FigureX.R_ & _SupFigureX.R_)

  It should be run in R with the models implemented in JAGS
