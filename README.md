# Penetrance
This project contains the core code used to estimate the penetrance of a gene.

## How it works
Since in many samples there were zero observations of CNVs in control cohorts, classical statistical methods using the binomial distribution or direct analysis of contingency tables cannot be applied and odds ratios are not defined. We therefore took a Bayesian approach, to derive a posterior distribution of likely values for the frequency of CNVs in case and control cohorts using the observed CNV frequencies from published data. We then simulated from this distribution to obtain empirical estimates of penetrance and their credible intervals.

We assumed a prior distribution for CNV frequency (p) for both cases and controls of p ~Beta(α, β), which is defined for 0<p<1, as required, with α, β > 0. Then, with n CNVs observed in N cases, and m CNVs observed in M controls (for n, m ≥ 0), the posterior distributions for CNV frequency in cases and controls are Beta(α+n, β+N-n) and Beta(α+m, β+M-m), respectively. We sampled pairs of CNV frequencies from these distributions, and calculate the penetrance, i.e. the probability of developing schizophrenia (disease, D) for individuals carrying the CNV (genotype, G), from

P(D/G) = P(G/D)P(D) / (P(G/D)P(D) + P(G/D')P(D'))

where D' represents controls, who do not
have schizophrenia, and P(D) the lifetime morbid risk for schizophrenia

## Usage
A web app to calculate the penetrance will be available soon.

## References
Penetrance for copy number variants associated with schizophrenia
Evangelos Vassos, David A. Collier, Simon Holden, Christine Patch, Dan Rujescu,David St Clair and Cathryn M. Lewis, June 2010


