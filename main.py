"""
TODO:
1. Still need to change code to include 2.5% and 97.5% quantiles and take median)
2. Add the confidence intervals
3. Include hidden cases of disease
4. Create GUI
5. Port to R?
"""
from scipy.stats import beta
import matplotlib.pyplot as plt
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

# Sample Data: 1q21.1 del from Kirov's data.
x = 1
y = 1
n = 17
N = 7918
m = 11
M = 46502

def probabilities(a,b):
    candidates = []
    # Display the Beta distibution
    fig, ax = plt.subplots(1, 1)
    x = np.linspace(beta.ppf(0.01, a, b),beta.ppf(0.99, a, b), 100)
    ax.plot(x, beta.cdf(x, a, b),'r-', lw=5, alpha=0.6, label='beta pdf')
    plt.show()
    # This is implementation only for the 50% quantile
    data_points = (np.asarray(beta.cdf(x, a, b)))
    for i,j in enumerate((range(len(data_points)))):
        if(round(data_points[j],1) == 0.5):
            # print(data_points[j])
            candidates.append(data_points[j])
    # Look for the value closes to 0.5 in the cdf of the beta distribution
    closest_value = min(candidates, key=lambda x:abs(x-0.5))
    # print(closest_value)
    return (beta.ppf(closest_value,a,b))

# For cases
a = x + n
b = y + N - n
prob_cases = (probabilities(a,b))

# For controls
A = x + m
B = y + M - m
prob_controls = (probabilities(A,B))

# Baysian Formula
lifetime_risk = 0.0072
penetrance = (prob_cases*lifetime_risk) / ((prob_cases*lifetime_risk)+(prob_controls*(1-(lifetime_risk))))
print(penetrance) # Should print 0.061 for L=1
