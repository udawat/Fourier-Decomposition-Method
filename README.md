# Fourier-Decomposition-Method
## Python implementation of Fourier Decomposition Method (FDM)

FDM is a nonlinear and nonstationary signal representation, decomposition and analysis method. It uses the Fourier and filter theory based zero-phase filtering for the decomposition of a signal into a set of suitable bands with desired cutoff frequencies (e.g. divide complete bandwidth of a signal into a set of sub-bands of equal bandwidth or dyadic subbands).

Perform FDM on given data `X`. The `fdm()` function expects X as a matrix of column vectors (horizontally-stacked signals). X may also be a single column vector (with one signal only). Here, `fs` is the sampling frequency. It generates subbands according to the given cutoff frequencies defined as a vector in `fc`.

```FIBFs = fdm(X, fs, fc)```

The function expects the input data as a Matrix of column vectors (horizontally-stacked signals). Using `data_type` as `'rows'` forces the function to process the input data matrix as row vectors (vertically-stacked signals).
```
FIBFs = fdm(X, fs, fc, data_type='columns')
FIBFs = fdm(X, fs, fc, data_type='rows')
```

To sort the output data from high to low frequency subbands (default behaviour). Use `sort_fc='ascend'` to change the default behaviour.
```
FIBFs = fdm(X, fs, fc, sort_fc='descend')
FIBFs = fdm(X, fs, fc, sort_fc='ascend')
```

To change the type of filter used from any of the following - `dft`, `dct`. (dct is used by default).
```
FIBFs = fdm(X, fs, fc, filter_type='dct')
FIBFs = fdm(X, fs, fc, filter_type='dft')
FIBFs = fdm(X, fs, fc, filter_type='fir')  # Implementation pending
FIBFs = fdm(X, fs, fc, filter_type='iir')  # Implementation pending
```

To generate a tiled layout figure showing all the different subbands. Default behaviour is to not show the plots.
```
FIBFs = fdm(X, fs, fc, plot_subbands=True)
FIBFs = fdm(X, fs, fc, plot_subbands=False)
```

Author: Original code for MATLAB written by Dr. Pushpendra Singh, Assistant Professor, National Institute of Technology, Hamirpur

Updated: Python Implementation, 25-Apr-2021 - Abhimanyu Singh Udawat, National Institute of Technology, Hamirpur

Modified: 03-Aug-2022 - Abhimanyu Singh Udawat

License:

For MATLAB based implementation of FDM, click here: https://www.researchgate.net/publication/363615241_An_efficient_MATLAB_code_for_faster_FDM_implementations_using_DFT_and_DCT

Please cite the following paper(s) when using this function:
**[The Fourier decomposition method for nonlinear and non-stationary time series analysis](https://royalsocietypublishing.org/doi/abs/10.1098/rspa.2016.0871)**

BibTeX Citation:
```
@article{doi:10.1098/rspa.2016.0871,
author = {Singh, Pushpendra  and Joshi, Shiv Dutt  and Patney, Rakesh Kumar  and Saha, Kaushik },
title = {The Fourier decomposition method for nonlinear and non-stationary time series analysis},
journal = {Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences},
volume = {473},
number = {2199},
pages = {20160871},
year = {2017},
doi = {10.1098/rspa.2016.0871},
URL = {https://royalsocietypublishing.org/doi/abs/10.1098/rspa.2016.0871},
eprint = {https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.2016.0871},
abstract = {For many decades, there has been a general perception in the literature that Fourier methods are not suitable for the analysis of nonlinear and non-stationary data. In this paper, we propose a novel and adaptive Fourier decomposition method (FDM), based on the Fourier theory, and demonstrate its efficacy for the analysis of nonlinear and non-stationary time series. The proposed FDM decomposes any data into a small number of ‘Fourier intrinsic band functions’ (FIBFs). The FDM presents a generalized Fourier expansion with variable amplitudes and variable frequencies of a time series by the Fourier method itself. We propose an idea of zero-phase filter bank-based multivariate FDM (MFDM), for the analysis of multivariate nonlinear and non-stationary time series, using the FDM. We also present an algorithm to obtain cut-off frequencies for MFDM. The proposed MFDM generates a finite number of band-limited multivariate FIBFs (MFIBFs). The MFDM preserves some intrinsic physical properties of the multivariate data, such as scale alignment, trend and instantaneous frequency. The proposed methods provide a time–frequency–energy (TFE) distribution that reveals the intrinsic structure of a data. Numerical computations and simulations have been carried out and comparison is made with the empirical mode decomposition algorithms.}
}
```
