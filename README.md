Implemented the paper "Diversity Gain for MIMO Neyman Pearson Signal Detection" which was published in IEEE Transactions on Signal Processing, March 2011

Overall Flow:

- Calculate the matrix A (Eq. 64)
- First we create the matrix Sn (Eq. 63)
- Then create C (Eq. 62)
- Then we create r, where zeta and w are generated using the random function (Eq 61)
- Then create the test statistic T (Eq 12)
- Monte Carlo for 100,000 

PS: The threshold values were chosen arbitrarily and may require fine-tuning to achieve a Probability of False Alarm of 0.001, ensuring alignment with the figures provided in the paper. (I did not perform this tuning as it took hours to run on my laptop.) Please feel free to raise any concerns regarding the implementation.



