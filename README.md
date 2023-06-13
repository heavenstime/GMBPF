# GMBPF (Gaussian Mixture Band Pass Filter)

## Files

- gmbpfPara.m : Parameters of GMBPF
- gmbpfApprox.m : Calculate coefficients to approximate GMBPF for given parameter
- gmbpfDesign.m : Optimize parameters
- gmbpf.m : Apply filter to a signal

## For the first time

- Make a folder named '''Data''' in the holder where programs exist. 

## Design a GMBPF

1. Set parameters in gmbpfPara.m
2. Excute gmbpfDesign.m (It calls gmbpfPara.m and gmbpfApprox.m)

### Note

- Until now, only P = 16 and P = 27 have been investigated in detail.
- They correspond to -60 dB and -100 dB stop band attenuations, respectively.
- Therefore, for other Ps, we must change the algorithm to decide the scan interval. 

### To show results of designed filter.

- Set flagPlot = 1 in gmbpfPara.m, and we can see the filter, its DTFT, its phase, and its time delay by plots. Furthermore, they are saved as text file. 
- As for file names, please see the source file. 

## Apply the GMBPF to a signal

1. Store the input signal into the varialbe `inSig` before execute gmbpf.m.
1. Excute gmbpf.m and the results are stored in `outSig`.
(The length of `outSig` is longer thane `inSig` by 2K. 

