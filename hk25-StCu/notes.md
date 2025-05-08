
# Mesoscale structure of stratocumulus (hk25-StCu): pre-hackathon intro notes

**Coordination**: Jakub Nowak (jakub.nowak@mpimet.mpg.de)



### Morphology metrics

Nonexhaustive list of cloud morphology metrics I encounterd in literature.


#### [Wood and Hartmann 2006](https://doi.org/10.1175/JCLI3702.1): statistics of liquid water path

Probability distributions
   - nondimensional homogeneity parameter $(\left<LWP\right>/\sigma_{LWP})^2$
   - skewness and kurtosis

2D power spectrum
   - The image is first detrended by removing the best-fit plane and then windowed using a Welch window.
   - The spectrum is a function of total wavenumber $k^2=k_x^2+k_y^2$
   - Two characteristic lengthscales are derived from the spectrum:
      - $\lambda_1$, peak wavelength, can be interpreted as the approximate diameter of mesoscale convective cells
      - $\lambda_2$, wavelength where the spectrum deviates from power law characteristic for smallest scales, can be interpreted as the typical size of clouds inside the cells.


#### [Wood et al. 2008](https://doi.org/10.1029/2007JD009371): variance of band-pass filtered IR brightness temperature

This metric was designed to detect pockets of open cells which exhibit higher mesoscale variability than their surrounding.

- A discrete wavelet transform is used as a method of band-pass filtering.
- Select a wavelet function and order. In practice, the order corresponds to the specific range of scales.
- Apply a discrete wavelet transform to each row of a 2D matrix and replace the values by the transform coefficients of the selected order.
- Repeat the operation for the columns.
- Take the square of the wavelet coefficients and smooth such a field with a running median of a rather large window.


#### [Koren et al. 2024](https://doi.org/10.1029/2024GL108435): cloud vs void chord length distributions (LvL)

- Consider a Bernoulli random matrix with M pixels where each pixel has an equal probability p of being cloudy.
- Rearrange the 2D matrix into a concatenated 1D vector by glueing sequentially. Do this flattening in both vertical and horizontal directions.
- The normalized random cloud-chord length distribution is $p^{L-1}(1-p)$ and the normalized random void-chord length distribution is $(1-p)^{L-1}p$.
- The deviation of an observed distribution from randomness is measured with a goodness-of-fit score based on the Kolmogorov-Smirnov test.
- LvL score consists of two components, for cloud and voids, forming the 2D LvL space.
- LvL fails when cloud fraction is 0 or 1!

There is a python [function](tools/LvL.py) provided in the [Koren et al. 2024 \[software\]](https://doi.org/10.34933/b7f2cded-40d3-4be9-bdc6-31b2694ca49c)


#### [Bagioli and Tompkins 2023](https://doi.org/10.1175/JAS-D-23-0103.1): Iorg and Lorg

This metric may be insuitable for stratocumulus clouds which are much connected objects but this can be verified.

Iorg
- Segment the binary field into connected objects.
- For each object, compute the distance from its centroid to a nearest neighbour.
- Derive the cumulative density function of nearest neighbour distances NNCDF(r).
- For a random field NNCDF(r) is a Weibull distribution $1-\exp(\lambda\pi r^2)$ (where $\lambda$ is mean density of objects).
- Plot the observed NNCDF against the NNCDF of a random field.
- Integrate the area under the graph to obtain Iorg.
- Iorg = 0.5 indicates randomness, Iorg > 0.5 clustering, Iorg < 0.5 regularity.

[comment]: # ( Iorg = 0.5 does not guarantee random distribution of objects because it integrates over all scales. There may be a superposition of clustered and regular subdistributions which appear at different scales. )

Lorg
- Instead of NNCDF(r), derive the distribution of number of neighbours NN(r) closer than r.
- Compute a radius of a disk which would include the same NN in the case of a random field, i.e. $\lambda\pi L^2=NN(r)$
- Plot $L(r)/r_{max}$ vs $r/r_{max}$ where $r_{max}$ is the maximum distance between the objects in the considered domain.
- Integrate the area under the graph to obtain Lorg.

There is a python [function](tools/ILorg.py) provided in https://github.com/giobiagioli/organization_indices.



### Conceptual mechanism

[Comstock et al. 2005](https://doi.org/10.1175/JAS3567.1), [van Zanten and Stevens 2005](https://doi.org/10.1175/JAS3611.1), [Wood et al. 2011](https://doi.org/10.5194/acp-11-2341-2011)


### Healpix grids

```python
nside = 2**zoom
ncells = 12*nside**2
```

| zoom | nside | res. (km) | ncells  |
| ----:| -----:| ---------:| ----------:|
|    0 | 	1 |	6519.6 | 12 |
|    1 | 	2 |	3259.8 | 48 |
|    2 | 	4 |	1629.9 | 192 |
|    3 |   	8 | 	815.0 | 768 |
|    4 |	16 | 	407.5 | 3,072 |
|    5 |	32 | 	203.7 | 12,288 |
|    6 |	64 | 	101.9 | 49,152 |
|    7 |   128 |  	50.9 | 196,608 |
|    8 |   256 |  	25.5 | 786,432 |
|    9 |   512 |  	12.7 | 3,145,728 |
|   10 |  1024 |   	6.4 | 12,582,912 |
|   11 |  2048 |   	3.2 | 50,331,648 |
|   12 |  4096 |   	1.6 | 201,326,592 |
