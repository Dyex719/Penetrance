## ABSTRACT

This project aims to automate the calculation involved in estimating the penetrance of a gene. This project will help researchers obtain the required penetrance by inputting values without having to worry about the internal mechanism behind the estimate.
The estimation is done using Bayesian probabilistic methods and the simulations to obtain these estimates was performed in Python using the SciPy package, with 2.5, 50 and 97.5% quantiles extracted to obtain the median penetrance, and its corresponding 95% credible intervals.

## INTRODUCTION

Penetrance is the probability of developing a particular disease given a particular genotype (genotypic model) or allele (allelic model). It is usually spoken in terms of lifetime risk, meaning the probability that the individual will develop that disease during their lifespan.
 
Conventional methods for estimating penetrance:
A typical study would look at everyone who has been observed with the given genotype, and ask how many have disease, or how many have disease by a certain age. This however presents many problems which tend to inflate one's estimate of penetrance.
 
In order to improve upon the traditional methods, researchers decided to use population based methods to study penetrance. The first key insight here is that the occurrence of a disease due to a genetic variant is directly proportional to the frequency of that variant in the population. The past few decades have seen a drastic increase in the amount of data collected for these types of studies and have helped to improve the accuracy of these results.

## CALCULATIONS

This project focuses on the paper written by Vassos et al. [1]  which describes all the procedure for estimating the penetrance using population based probabilistic techniques and simulations.

Since in some samples there were zero observations of CNVs in control cohorts, classical statistical methods using the binomial distribution or direct analysis of contingency tables cannot be applied and odds ratios are not defined.

 The researchers therefore took a Bayesian approach, to derive a posterior distribution of likely values for the frequency of CNVs in case and control cohorts using the observed CNV frequencies from published data. Simulations were then performed from this distribution to obtain empirical estimates of penetrance and their credible intervals.

A prior distribution for CNV frequency (p) for both cases and controls of p ~Beta(α, β), which is defined for 0<p<1, as required, with α, β > 0. Then, with n CNVs observed in N cases, and m CNVs observed in M controls (for n, m ≥ 0), the posterior distributions for CNV frequency in cases and controls are Beta(α+n, β+N-n) and Beta(α+m, β+M-m), respectively. Pairs of CNV frequencies were sampled from these distributions, and were used calculate the penetrance, i.e. the probability of developing schizophrenia (disease, D) for individuals carrying the CNV (genotype, G), from

![Results](/images/Bayes_Formula.jpg)

where D' represents controls, who do not have schizophrenia, and P(D) the lifetime morbid risk for schizophrenia

Method to calculate the penetrance as described in [1]:

At first, the above formula replacing P(G| D) with n/N (CNV frequency in
cases) and P(G | D) with m/M (CNV frequency in controls). This approach assumes that
the controls of the original studies did not include cases with schizophrenia.

 However, since not all the controls were screened for schizophrenia and some controls may develop the disease after their recruitment, the possibility that the control groups may include “hidden” cases with schizophrenia was examined.

 Let L be the proportion of controls screened absent for schizophrenia, and assume, conservatively, that the remaining proportion of controls (1-L) carry the general population risk for schizophrenia. Thus, the total number of cases was estimated to be N+M*(1-L)*P(D), the number of screened controls M- M*(1-L)*P(D), the expected number of CNVs in cases n+M*(1-L)*P(D)*n/N and the number of CNVs in controls m- M*(1-L)*P(D)*n/N.

 Using these formulae the penetrance was estimated for various values of L including the extreme values of L=1 (no schizophrenia cases hidden in controls, as assumed previously) and L=0 (general population controls, unscreened for schizophrenia).

For the estimation of penetrance of schizophrenia, the conservative median global value for lifetime morbid risk (P(D) in the formula mentioned above) of 0.72% was adopted from a comprehensive meta-analysis. Simulation was performed in Python using the SciPy package, with 2.5, 50 and 97.5% quantiles extracted to obtain the median penetrance, and its 95% credible intervals.

Method to calculate the corresponding credible intervals is described in [2]:

The upper bound of the 95% credible interval is calculated based on the upper 95%CI of the binomial distribution for case allele frequency and the lower 95%CI for controls. Conversely, the lower bound of penetrance is based on the lower bound of case allele frequency and the upper bound of control allele frequency.

## WEB APPLICATION

The webapp was created using Flask and PythonAnywhere.
 
Flask is a micro web framework written in Python and based on the Werkzeug toolkit and Jinja2 template engine.
 
PythonAnywhere is an online integrated development environment (IDE) and web hosting service based on the Python programming language.

Web application available at: http://calpen.pythonanywhere.com/


## RESULTS

The following shows the results using the data from paper[1]:
![Results](/images/Results1.jpg)
 
The following shows the results using the data from paper[3]:
![Results](/images/Results2.jpg)


## REFERENCES


1.	Penetrance for copy number variants associated with schizophrenia
Evangelos Vassos, David A. Collier, Simon Holden, Christine Patch, Dan Rujescu, David St Clair and Cathryn M. Lewis.

2.	The Penetrance of Copy Number Variations for Schizophrenia and Developmental Delay
George Kirov, Elliott Rees, James T.R. Walters, Valentina Escott-Price, Lyudmila Georgieva, Alexander L. Richards, Kimberly D. Chambert, Gerwyn Davies, Sophie E. Legge, Jennifer L. Moran, Steven A. McCarroll, Michael C. O’Donovan, and Michael J. Owen.

3.	Estimates of penetrance for recurrent pathogenic copy-number variations
Jill A. Rosenfeld, Bradley P. Coe, Evan E. Eichler, Howard Cuckle and Lisa G. Shaffer.





