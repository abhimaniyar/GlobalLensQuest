#  GlobalLensQuest

Comparison of various CMB lensing quadratic estimators in terms of their noise power spectrum for any CMB experiment:

* Global minimum variance quadratic estimator based on an alternate derivation
* Standard Hu & Okamoto 2002 quadratic estimators from temperature and polarization. This is not the most-optimum quadratic estimator contrary to what was thought before
* Okamoto and Hu 2003 full-sky version of Hu and Okamoto estimators. Has a provision of setting C^{TE}_\ell = 0 to have a separation in the configuration space
* Sub-optimal quadratic estimator. This was used in the Planck 2018 lensing analysis where they set C^{TE}_\ell = 0 in their covariance matrix. Results in an estimator suboptimal even than Hu and Okamoto 2002

Requires numpy for computation and matplotlib for plotting

You can just clone the repository, and do all the calculations as in driver.py:
```
python driver.py
```
This code uses some of the structure from Manu Schaan's code ForQuE found at https://github.com/EmmanuelSchaan/ForQuE
If you have any comments, suggestions or questions, please do not hesitate to contact me: abhishek.maniyar@nyu.edu

