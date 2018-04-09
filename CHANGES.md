Version 0.9.0 (2018-04-09):
--------------------------
  - Incorporate a coverage likelihood that corrects for the PHMM emission
    likelihoods preference. This effectively acts as a regularizer that
    incorporates more information about how we expect reads to be distributed
    across a diploid.

Version 0.8.0 (2018-01-08):
--------------------------
 - First public release. I have arbitrarily chosen version 0.8 to reflect the
   maturity of the project. In the current state I believe that the tool will
   give valid HLA type inference in most (90%) cases, and there are very few
   bugs in the implementation and logic described in the paper.
 - Future releases, before version 1.0, will address other parameters that we
   believe can help with robust typing such as a uniform/consistent coverage
   between chosen diploid pairs.
