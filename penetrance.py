from scipy.stats import beta
import scipy.stats as st
import numpy as np
from scipy import stats, optimize
from sys import float_info

def check_input(prompt, default):
    while True:
        try:
            value = int(input(prompt) or default)
        except ValueError:
            print("Please enter a positive natural number.\n")
            continue
        if value < 0:
            print("Sorry, your response must not be negative.\n")
            continue
        else:
            break
    return value

def check_input_float(prompt, default):
    while True:
        try:
            value = (float(input(prompt) or default))/100
        except ValueError:
            print("Please enter a number in the form of the percentage affected.\n")
            continue
        if value < 0:
            print("Sorry, your response must not be negative.\n")
            continue
        else:
            break
    return value

def Bayesian(prob_cases,prob_controls, hidden = 0):
    lifetime_risk = 0.0072
    return (prob_cases*lifetime_risk) / ((prob_cases*lifetime_risk)+(prob_controls*(1-(lifetime_risk))))

def probabilities(a,b,quantiles):
    candidates = []
    # Define the Beta distibution
    # Increase this value to increase the precision
    x = np.linspace(beta.ppf(0.001, a, b),beta.ppf(0.999, a, b), 10000)
    data_points = (np.asarray(beta.cdf(x, a, b)))

    for i,j in enumerate((range(len(data_points)))):
        if(round(data_points[j],3) == quantiles):
            # print(data_points[j])
            candidates.append(data_points[j])

    # Look for the value closest to 0.5 in the cdf of the beta distribution
    closest_value = min(candidates, key=lambda x:abs(x-quantiles))
    return (beta.ppf(closest_value,a,b))

def proportion_confint(count, nobs, alpha=0.05, method='normal'):
    q_ = count * 1. / nobs
    alpha_2 = 0.5 * alpha

    if method == 'wilson':
        crit = stats.norm.isf(alpha / 2.)
        crit2 = crit**2
        denom = 1 + crit2 / nobs
        center = (q_ + crit2 / (2 * nobs)) / denom
        dist = crit * np.sqrt(q_ * (1. - q_) / nobs + crit2 / (4. * nobs**2))
        dist /= denom
        ci_low = center - dist
        ci_upp = center + dist
    else:
        raise NotImplementedError('method "%s" is not available' % method)
    return (ci_low, ci_upp)

def penetrance_confint(ac_case, n_case, ac_control, n_control, baseline_risk):
    case_confint =  proportion_confint(count=ac_case,nobs=2*n_case,method='wilson')
    control_confint =  proportion_confint(count=ac_control,nobs=2*n_control,method='wilson')

    lower_bound = Bayesian(case_confint[0],control_confint[1])
    upper_bound = Bayesian(case_confint[1],control_confint[0])

    return (lower_bound,upper_bound)

if __name__ == "__main__":
    x = check_input("Enter beta dist parameter a: (Default = 1)\n",1)
    y = check_input("Enter beta dist parameter b: (Default = 1)\n",1)
    n = check_input("Enter the number of CNVs observed in cases: (Default = 8)\n",8)
    N = check_input("Enter the total number of patients in cases: (Default = 5089)\n",5089)
    m = check_input("Enter the number of CNVs observed in controls: (Default = 6)\n",6)
    M = check_input("Enter the total number of patients in controls: (Default = 38884)\n",38884)
    P = check_input_float("Enter the proportion of controls screened for the disease (Default = 0.72%):\n",0.72)
    L = 1

    n = n + M*(1-L)*P*n/N
    N = N + M*(1-L)*P
    m = max(0, m - M*(1-L)*P*n/N)
    M =  M - M*(1-L)*P

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
    lower,upper = penetrance_confint(n,N,m,M,P)

    print("The penetrance is %.3f\n" % median_penetrance) # Should print 0.067 for L=1
    print("The confidence interval ranges from: %.3f to %.3f\n" % (lower,upper))

main.py
Displaying main.py.
