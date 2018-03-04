"""
TODO:
1. Still need to change code to include 2.5% and 97.5% quantiles and take median DONE!
2. Add the confidence intervals DONE!
3. Include hidden cases of disease DONE!
4. Create web app (ongoing)
"""
from scipy.stats import beta
import scipy.stats as st
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np

# print("Enter beta dist parameter a:")
# x = int(input())
# print("Enter beta dist parameter b):")
# y = int(input())
# print("Enter the number of CNVs observed in cases:")
# n = int(input())
# print("Enter the total number of patients in cases:")
# N = int(input())
# print("Enter the number of CNVs observed in controls:")
# m = int(input())
# print("Enter the total number of patients in controls:")
# M = int(input())
# print("Enter the proportion of controls screened for schizophrenia:")
# M = int(input())

# Sample Data: 22q11 del from Kirov's data.
x = 1
y = 1

# Includes the hidden case factor of unscreened controls developing schizophrenia.
n = 8
N = 5089
m = 6
M = 38884
P = 0.0072
L = 0.5

n = n + M*(1-L)*P*n/N
N = N + M*(1-L)*P
m = max(0, m - M*(1-L)*P*n/N)
M =  M - M*(1-L)*P

def mean_confidence_interval(data, confidence=0.95):
    return st.t.interval(confidence, len(data)-1, loc=np.mean(data), scale=st.sem(data))

def Bayesian(prob_cases,prob_controls, hidden = 0):
    lifetime_risk = 0.0072
    return (prob_cases*lifetime_risk) / ((prob_cases*lifetime_risk)+(prob_controls*(1-(lifetime_risk))))

def plot(x):
    fig, ax = plt.subplots(1, 1)
    ax.plot(x, beta.cdf(x, a, b),'r-', lw=5, alpha=0.6, label='beta pdf')
    plt.show()

def probabilities(a,b,quantiles):
    candidates = []
    # Define the Beta distibution
    # Increase this value to increase the precision
    x = np.linspace(beta.ppf(0.001, a, b),beta.ppf(0.999, a, b), 10000)

    data_points = (np.asarray(beta.cdf(x, a, b)))
    # credible_int = bayes_mvs(data_points,0.95)
    # credible_int = mean_confidence_interval(x,0.95)
    # print(credible_int)
    # quit()
    for i,j in enumerate((range(len(data_points)))):
        if(round(data_points[j],3) == quantiles):
            # print(data_points[j])
            candidates.append(data_points[j])

    # Look for the value closest to 0.5 in the cdf of the beta distribution
    closest_value = min(candidates, key=lambda x:abs(x-quantiles))
    # print(closest_value)
    return (beta.ppf(closest_value,a,b))

# For cases
a = x + n
b = y + N - n

# For controls
A = x + m
B = y + M - m

# Baysian Formula
penetrance_5 = Bayesian(probabilities(a,b,0.500),probabilities(A,B,0.500))
penetrance_25 = Bayesian(probabilities(a,b,0.025),probabilities(A,B,0.025))
penetrance_975 = Bayesian(probabilities(a,b,0.975),probabilities(A,B,0.975))
median_penetrance = np.median([penetrance_5,penetrance_25,penetrance_975])

print("%.3f" % median_penetrance) # Should print 0.553 for L=1
