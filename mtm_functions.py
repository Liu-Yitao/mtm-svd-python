# Script for MultiTaper Method-Singular Value Decomposition (MTM-SVD) with Monte Carlo test in python

# ------------------------------------------------------------------

# This script is a modified version of the Python function developed by
# Mathilde Jutras at McGill University, Canada[1]. 
# You can find the original Python code here: 
# https://github.com/mathildejutras/mtm-svd-python

# This script was adapted by Yitao Liu at Nanjing University of Information Science & Technology, China
# Copyright (C) 2021, Yitao Liu
# and is available under the GNU General Public License v3.0

# The script may be used, copied, or redistributed as long as it is cited as follow:
# Liu-Yitao, & Mathilde Jutras. (2021). Liu-Yitao/mtm-svd-python: MTM-SVD with Monte Carlo test (v1.1.0). Zenodo. https://doi.org/10.5281/zenodo.5774584

# This software may be used, copied, or redistributed as long as it is not 
# sold and that this copyright notice is reproduced on each copy made. 
# This routine is provided as is without any express or implied warranties.

# Questions or comments to:
# Yitao Liu, liuyitao97@outlook.com

# Last update:
# Dec 2021

# ------------------------------------------------------------------

# The script is structured as follows:

# In the main script is found in mtm-svd-python.py
# In the first section, the user can load the data,
# assuming the outputs are stored in a netcdf format.
# In the secton section, functions are called to calculate the spectrum
# The user will then be asked for which frequencies he wants to plot 
# the spatial patterns associated with the variability.
# In the third section, the spatial patterns are plotted and saved

# The required functions are found in mtm_functions.py

# ------------------------------------------------------------------

# Python Package needed:
# - numpy
# - scipy
# - xarray (read the netcdf file)
# - matplotlib (not necessary, just for plotting)

# You can install the needed Python packages by conda,with the command below
# ```
# conda install -c conda-forge numpy scipy xarray matplotlib
# ```

# [1] Mathilde Jutras. (2020, July 6). mathildejutras/mtm-svd-python: v1.0.0-alpha (Version v1.0.0). Zenodo. http://doi.org/10.5281/zenodo.3932319

# ------------------------------------------------------------------

import numpy as np
from numpy.random import shuffle

def mtm_svd_lfv(ts2d, psi, kk, dt):
    """
    Calculate the local fractional variance spectrum (LFV) using the MultiTaper Method-Singular Value Decomposition (MTM-SVD).

    Parameters:
    - ts2d: 2D array of shape (p, n), where p is the number of grid points and n is the number of time steps.
    - psi: Discrete Prolate Spheroidal Sequences
    - kk: Number of Slepian tapers to use.
    - dt: Time step size.

    Returns:
    - fr: Array of frequencies.
    - lfvs: Array of LFV values corresponding to each frequency.
    """

    # Compute spectrum at each grid point
    p, n = ts2d.shape

    # Remove the mean and divide by std
    vm = np.nanmean(ts2d, axis=0)  # Calculate the mean along the time axis
    ts2d = ts2d - vm  # Subtract the mean from ts2d
    vs = np.nanstd(ts2d, axis=0)  # Calculate the standard deviation along the time axis
    ts2d = np.divide(ts2d, vs, where=vs!=0)  # Divide ts2d by the standard deviation, avoid division by zero
    ts2d = np.nan_to_num(ts2d)  # Replace NaN values with 0

    npad = 2**int(np.ceil(np.log2(abs(p))) + 2)  # Calculate the padding size for the FFT
    nf = int(npad / 2)  # Calculate the number of frequencies
    ddf = 1./(npad*dt)  # Calculate the frequency resolution
    fr = np.arange(0, nf) * ddf  # Generate the array of frequencies

    # Get the matrix of spectrums
    psimats = np.array([np.multiply(psi[k, :, np.newaxis], ts2d) for k in range(kk)])  # Use broadcasting instead of repmat
    nev = np.fft.fft(psimats, n=npad, axis=1)  # Perform the FFT along the time axis
    nev = np.fft.fftshift(nev, axes=(1))  # Shift the FFT output to center the frequencies
    nev = nev[:, nf:, :]  # Remove the negative frequencies

    # Calculate svd for each frequency
    lfvs = np.zeros(nf) * np.nan  # Initialize an array to store the LFV values
    for j in range(nf):
        U, S, V = np.linalg.svd(nev[:, j, :], full_matrices=False)  # Perform the SVD on the matrix at each frequency
        lfvs[j] = S[0]**2 / (np.nansum(S[1:])**2)  # Calculate the LFV value for the current frequency

    return fr, lfvs

# Function 2) Calculate the confidence interval of the LFV calculations
def monte_carlo_test(index_with_space, niter, sl, len_freq, psi, kk, dt):
	"""
	Perform Monte Carlo test to calculate the confidence interval of the LFV calculations.

	Parameters:
	- index_with_space: 2D array of shape (p, n), where p is the number of grid points and n is the number of time steps.
	- niter: int, Number of iterations for the Monte Carlo test.
	- sl: list, Significance level for the confidence interval.
	- len_freq: Number of frequencies.
	- psi: Discrete Prolate Spheroidal Sequences
	- kk: Number of Slepian tapers to use.
	- dt: Time step size.

	Returns:
	- freq_rd: Array of frequencies.
	- conflev: Array of LFV values representing the confidence interval.

	"""

	# create file to store all random lfv
	lfv_mc = np.zeros((niter, len_freq))

	# calculate all random lfc and store in lfv_mc[num_of_mc, num_of_freq]
	for ii in range(niter):
		if (ii % 100) == 0:
			print(f'niter = {ii}')
		shuffle(index_with_space)  # randomly shuffle the index_with_space array
		_, lfv_rd = mtm_svd_lfv(index_with_space, psi, kk, dt)  # calculate the LFV using the shuffled index_with_space
		lfv_mc[ii, :] = lfv_rd  # store the LFV values in the lfv_mc array

	lfv_mc_sort = np.sort(lfv_mc, axis=0)  # sort the LFV values in ascending order

	num = np.rint((1 - np.asarray(sl)) * niter).astype(np.int64)  # calculate the number of LFV values to include in the confidence interval

	conflev = lfv_mc_sort[num, :]  # select the LFV values corresponding to the confidence interval

	return conflev


# @jit
# @jit
def envel(ff0, iif, fr, dt, ddf, p, kk, psi, V):
	"""
	Calculate the envelope of a signal at a specific frequency using the MTM-SVD method.

	Parameters:
	- ff0: The target frequency.
	- iif: The index of the target frequency in the frequency array.
	- fr: Array of frequencies.
	- dt: Time step size.
	- ddf: Frequency resolution.
	- p: Number of grid points.
	- kk: Number of Slepian tapers.
	- psi: Array of Slepian tapers.
	- V: Matrix of singular vectors.

	Returns:
	- env: Array of envelope values.

	"""

	ex = np.ones(p)  # Initialize an array of ones with length p
	df1 = 0  # Initialize the frequency difference

	c0 = 1
	s0 = 0

	c = [c0]
	s = [s0]
	cs = np.cos(2. * np.pi * df1 * dt)
	sn = np.sin(2. * np.pi * df1 * dt)
	for i in range(1, p):
		c.append(c[i - 1] * cs - s[i - 1] * sn)  # Calculate the cosine term
		s.append(c[i - 1] * sn + s[i - 1] * cs)  # Calculate the sine term

	cl = np.ones(p)  # Initialize an array of ones with length p
	sl = np.zeros(p)  # Initialize an array of zeros with length p

	d = V[0, :]  # Get the first row of the singular vectors
	d = np.conj(d) * 2  # Take the complex conjugate and multiply by 2
	if iif == 1:
		d = V[0, :]  # Get the first row of the singular vectors
		d = np.conj(d)  # Take the complex conjugate

	g = []
	for i0 in range(kk):
		cn = [complex(psi[i0, i] * c[i], -psi[i0, i] * s[i]) for i in range(len(s))]  # Calculate the complex number
		g.append(ex * cn)  # Multiply the complex number by ex and append to the list
	g = np.array(g).T  # Convert the list to a numpy array and transpose

	za = np.conj(sum(g))  # Calculate the complex conjugate of the sum of g

	[g1, qrsave1] = np.linalg.qr(g)  # Perform QR decomposition on g
	dum1 = np.linalg.lstsq(np.conj(qrsave1), np.linalg.lstsq(np.conj(qrsave1.T), d)[0])[0].T  # Solve the linear equation
	amp0 = sum(np.conj(za) * dum1)  # Calculate the amplitude
	dum2 = np.linalg.lstsq(np.conj(qrsave1), np.linalg.lstsq(np.conj(qrsave1.T), za)[0])[0].T  # Solve the linear equation
	amp1 = sum(np.conj(za) * dum2)  # Calculate the amplitude
	amp0 = amp0 / amp1  # Normalize the amplitude
	sum1 = sum(abs(d) ** 2)  # Calculate the sum of squared absolute values of d
	d = d - za * amp0  # Subtract the product of za and amp0 from d
	sum2 = sum(abs(d) ** 2)  # Calculate the sum of squared absolute values of d
	env0 = np.linalg.lstsq(np.conj((qrsave1.T)), d.T)[0].T  # Solve the linear equation
	env = np.matmul(g1, env0.T)  # Multiply g1 and env0 and transpose
	env = env + amp0 * np.ones(len(c))  # Add amp0 to env

	return env


# Function 3) Reconstruct the spatial patterns associated with peaks in the spectrum
# @jit

def mtm_svd_recon(ts2d, psi, kk, dt, fo):
    """
    Perform multitaper singular value decomposition (SVD) reconstruction.
    Args:
        ts2d (ndarray): 2D array of time series data.
        psi (ndarray): Discrete Prolate Spheroidal Sequences.
        kk (int): Number of tapers to use.
        dt (float): Time step between samples.
        fo (list): List of frequencies of interest.
    Returns:
        tuple: A tuple containing the following:
            - vexp (list): List of variance explained for each frequency.
            - totvarexp (list): List of total variance explained for each frequency.
            - iis (list): List of indices corresponding to the closest frequency for each frequency of interest.
    """
    imode = 0  # Initialize the mode index

    # Compute spectrum at each grid point
    p, n = ts2d.shape  # Get the shape of the input time series data

    # Remove the mean and divide by std
    vm = np.nanmean(ts2d, axis=0)  # Calculate the mean along the time axis
    ts2d = ts2d - vm  # Subtract the mean from ts2d
    vs = np.nanstd(ts2d, axis=0)  # Calculate the standard deviation along the time axis
    ts2d = np.divide(ts2d, vs, where=vs != 0)  # Divide ts2d by the standard deviation, avoid division by zero
    ts2d = np.nan_to_num(ts2d)  # Replace NaN values with 0

    npad = 2**int(np.ceil(np.log2(abs(p))) + 2)  # Calculate the padding size for the FFT
    nf = int(npad / 2)  # Calculate the number of frequencies
    ddf = 1. / (npad * dt)  # Calculate the frequency resolution
    fr = np.arange(0, nf) * ddf  # Generate the array of frequencies

    # Get the matrix of spectrums
    psimats = np.array([np.multiply(psi[k, :, np.newaxis], ts2d) for k in range(kk)])  # Use broadcasting instead of repmat
    nev = np.fft.fft(psimats, n=npad, axis=1)  # Perform the FFT along the time axis
    nev = np.fft.fftshift(nev, axes=(1))  # Shift the FFT output to center the frequencies
    nev = nev[:, nf:, :]  # Remove the negative frequencies

    # Initialize output matrices
    S = np.ones((kk, len(fo))) * np.nan  # Initialize the matrix S with NaN values
    vexp = []  # Initialize an empty list to store the variance explained for each frequency
    totvarexp = []  # Initialize an empty list to store the total variance explained for each frequency
    iis = []  # Initialize an empty list to store the indices corresponding to the closest frequency for each frequency of interest

    D = vs  # Assign the standard deviation array to D

    for i1 in range(len(fo)):
        # closest frequency
        iif = (np.abs(fr - fo[i1])).argmin()  # Find the index of the closest frequency to the current frequency of interest
        iis.append(iif)  # Append the index to the list
        ff0 = fr[iif]  # Get the closest frequency
        print('( %i ) %.2f cyclesyr | %.2f yr' % (iif, ff0, 1 / ff0))  # Print the index, frequency in cycles per year, and frequency in years

        U, S0, Vh = np.linalg.svd(nev[:, iif, :].T, full_matrices=False)  # Perform the SVD on the matrix at the closest frequency
        V = Vh  # Transpose the Vh matrix
        S[:, i1] = S0  # Store the singular values in the S matrix

        env1 = envel(ff0, iif, fr, dt, ddf, p, kk, psi, V)  # Calculate the envelope using the envel function

        # Calculate cosine and sine terms
        t = np.arange(p)
        c = np.cos(2 * np.pi * ff0 * dt * t)
        s = np.sin(2 * np.pi * ff0 * dt * t)
        CS = c + 1j * s
        CS = np.conj(CS)  # Take the complex conjugate of the complex number

        # Reconstructions
        R = np.real(D * S[imode, i1] * np.outer(U[:, imode], CS * env1).T)  # Perform the reconstruction

        vsr = np.var(R, axis=0)  # Calculate the variance along the time axis
        vexp.append(vsr / (vs ** 2) * 100)  # Calculate the variance explained and append to the list
        totvarexp.append(np.nansum(vsr) / np.nansum(vs ** 2) * 100)  # Calculate the total variance explained and append to the list

    return vexp, totvarexp, iis  # Return the variance explained, total variance explained, and indices
