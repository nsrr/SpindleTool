# SpindleTool
GUI for spindle detection from data in NSRR format using wavelets or bandpass filtering.

Adapted from methods described in Purcell, S.M., Manoach, D.S., Demanuele, C., Cade, B.E., Mariani, S., Cox, R., Panagiotaropoulou, G., Saxena, R., Pan, J.Q., Smoller, J.W. and Redline, S., 2017. Characterizing sleep spindles in 11,630 individuals from the National Sleep Research Resource. Nature Communications, 8.

Your data will need to be in the same format as NSRR data: EDF for signals and XML for annotations

The two alternative methods used for spindle detection are:
- wavelet-based --> decribed in E. J. Wamsley, M. A. Tucker, A. K. Shinn, K. E. Ono, S. K. McKinley, A. V. Ely, D. C. Goff, R. Stickgold, and D. S. Manoach, “Reduced Sleep Spindles and Spindle Coherence in Schizophrenia: Mechanisms of Impaired Memory Consolidation?,” Biol. Psychiatry, vol. 71, no.2, pp. 154–161, Jan. 2012.
A Morlet wavelet is used for a first detection, combined with a threshold criterion. A second threshold is used for the calculation of the extended portion of the spindles.

- bandpass-based --> described in M. D. Ferrarelli Fabio, P. D. Huber Reto, M. D. Peterson Ph. D.,Michael, M. D. Massimini Ph. D.,Marcello, B. S. Murphy Michael, B. S. Riedner Brady, A. Watson, M. D. Bria Pietro, and M. D. Tononi Ph. D.,Giulio, “Reduced Sleep Spindle Activity in Schizophrenia Patients,” Am. J. Psychiatry, vol. 164, no. 3, pp. 483–492, Mar. 2007.
A bandpass fir filter is used to obtain power in the sigma band. A threshold is superimposed to the results.

Please report bugs to sara.mariani7@gmail.com
