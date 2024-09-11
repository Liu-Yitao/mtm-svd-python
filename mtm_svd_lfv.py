import numpy as np
from scipy.linalg import svd
from scipy.signal.windows import dpss

def mtm_svd_lfv(sh, nw, k, dt, vw, npad):
    """
    Determine the local fractional variance spectrum LFV.

    Parameters:
        sh : ndarray
            Time series matrix [rows=time columns=space]
        nw : int
            Time-frequency bandwidth parameter
        k : int
            Number of Slepian Tapers.
        dt : float or ndarray
            Time interval (cycles/year) or date vector.
        vw : ndarray
            Weights on dataseries.
        npad : int
            Zero-padding factor.

    Returns:
        LFV : dict
            Dictionary containing the LFV spectrum and related parameters.
    """
    print('MTM-SVD LOCAL FRACTIONAL VARIANCE SPECTRUM')
    print('Preprocessing...')
    
    n, p = sh.shape

    # # Time interval
    # if len(dt) > 1:
    #     lt = (dt[-1] - dt[0]) / 365.25  # Length of the series in years
    #     dt = lt / len(dt)

    # Remove mean and normalize by standard deviation
    vm = np.nanmean(sh, axis=0)
    sh -= vm
    vs = np.nanstd(sh, axis=0)
    vs[vs == 0] = 1  # Avoid division by zero
    sh /= vs

    # Initialize time series with uniform weights. Can use other weights (e.g., area).
    if vw == 0:
        vw = np.ones(p)

    # Weighting by the weights
    sh *= vw

    # Fill missing data with linear interpolation
    for i in range(p):
        idx = np.isnan(sh[:, i])
        if np.any(idx):
            ai = np.flatnonzero(~idx)
            if len(ai) > 1:
                sh[idx, i] = np.interp(np.flatnonzero(idx), ai, sh[ai, i])
            sh[np.isnan(sh)] = 0

    # Slepian tapers
    psi = dpss(n, nw, k)

    # Padding with zeros
    if npad == 0:
        npad = 2 ** (int(np.log2(n)) + 1)

    # Display spectral estimation parameters
    print(f'Slepian Tapers : {k}')
    print(f'Bandwidth parameter : {nw}')
    print(f'Padding : {npad}')
    print('Performing...')

    nf = npad // 2  # Half spectrum dimensions
    ddf = 1 / (npad * dt)  # Sampling frequency in npad intervals
    fr = (np.arange(nf) * ddf)  # Frequency vector

    # Perform spectral estimation
    nev = []
    for i1 in range(k):
        psimat = np.repeat(psi[:, i1:i1+1], p, axis=1) * sh
        psimat = np.fft.fft(psimat, npad, axis=0)[:nf]
        skn = f'sk{i1}'
        nev.append(f'{skn}(j1, :).T')
        print(f'> Spectral estimated {i1}: ')

    # SVD of the spectrum matrix for different frequencies
    LFVS = np.full(nf, np.nan)
    for j1 in range(nf):
        M = np.array([eval(nev[i1]) for i1 in range(k)])
        S = svd(M, full_matrices=False)[1]
        LFVS[j1] = S[0] ** 2 / sum(S[1:] ** 2)

    vfbw = 2 * nw * (1 / dt) * (1 / n)

    # Store results in a structured dictionary
    LFV = {
        'name': 'Local Fractional Variance Spectrum',
        'timeinterval': dt,
        'bandwidth': nw,
        'tapers': k,
        'padding': npad,
        'varfreqbandwidth': vfbw,
        'spectrum': LFVS,
        'specdomain': fr
    }

    print(f'> Local Fractional Variances : ')
    print('Done.')

    return LFV


# Example usage:

sh = np.random.rand(100, 5)  # Replace with actual time series data
LFV = mtm_svd_lfv(sh, 2, 3, 0.09, 0, 512)