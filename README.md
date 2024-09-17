# TimeSeries_BA
**Start date:** 2024.04.12  
**Archived date:** 2024.06.07

This is an archive of a pilot experiment of DNAm&TSA

## Introduction

In the field of aging research, the concept that an organism's degree of aging is exactly equal to the number of years they have lived, that is, their age, has been refuted. For example, in a population, people of the same age vary in their degree of aging (in fact, it is inaccurate to use the term 'aging' in such a broad sense—what does 'aging' really mean?). Even within the same person, different organs age at different rates. A more extreme example is that people with progeria age far more rapidly than their chronological age would suggest.Researchers have developed various algorithms to measure a general concept of aging using the concept of biological age. One major class of algorithms is based on combinations of blood test biomarkers, and another major class is based on DNA methylation data (different parts of the DNA strand have varying degrees of methylation, which we measure to determine the methylation probability at each site, then use data from numerous sites to train models).

This paper treats a person's biological age (calculated using DNA methylation methods, including the Horvath1, Hannum2, Lin3, and Ying_caus4 clocks) as a simulated time series. The starting point of the biological age-based pseudo time series is t=1, with one-year intervals, meaning that a person's biological age of 1, 2, 3, 4, 5, ..., 120, ... represents t=1, 2, 3, ... 120. The series does not start at 0 because, in the ground zero model, it is believed that the initial embryo and the biological age at birth are not 0 but follow a smooth, approximate V-shaped curve, reaching the lowest point in the mid-embryonic stage. The methylation β value of a given site is treated as a probability. Using a 450k methylation array as an example, with over 450,000 observation sites, considering all sites would result in panel data. The significance of this modeling might be to study the stochastic mechanisms that govern changes in human DNA methylation and whether an individual's DNA methylation profile is time-dependent.

One practical issue is that cross-sectional DNA methylation data from large populations is relatively easy to obtain, whereas long-term longitudinal data is not. Since this project attempts TSA learning, I have made significant simplifications by averaging the data of individuals of the same biological age (e.g., 20 years old). This allows n individuals to represent an "average 20-year-old," and this approach is extrapolated to each year to construct a statistically "average person" methylation panel from age 1 to 120. Given that the central limit theorem applies here, this should approximate the methylation levels of a specific individual at each age. As a result, the methodology of this project is somewhat closer to an ecological study, which inherently comes with the common statistical pitfall of ecological research—limited extrapolability.

This paper utilizes ten datasets from the Gene Expression Omnibus: GSE41169, GSE40279, GSE51057, GSE42861, GSE51032, GSE73103, GSE69270, GSE64495, GSE30870, and GSE52588, totaling 3392 samples, all of which are Homo sapiens samples. The detection method used is the Illumina HumanMethylation450 BeadChip, with standardized 485,577 CpG sites, and the data have passed heterogeneity tests on the Gene Expression Omnibus platform.

## Methodology

1. **Establishing the time series using aging clocks**
2. **Descriptive statistics, handling of missing values, and outliers**
3. **Testing the stationarity of the original series**; if the original series is non-stationary, use differencing to remove the trend
4. **Determining if the series exhibits seasonality (periodicity)**; if there is seasonal (periodic) fluctuation, use seasonal (periodic) differencing to remove seasonality (periodicity)
5. **Identifying the model**
6. **Determining the order of the model**
7. **Estimating the parameters of the model**
8. **Testing the goodness of fit of the model** by conducting a white noise test on the residual series to determine if it is a white noise series
9. **Providing the model's prediction results and plotting the trend forecast graph**


