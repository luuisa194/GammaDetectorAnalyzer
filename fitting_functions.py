import numpy as np

# Constant for 2π
TWO_PI = 2 * np.pi

def normalized_gaussian(x, mu, sigma, amplitude):
    """
    Computes a normalized Gaussian function.
    """
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2) / np.sqrt(TWO_PI * sigma ** 2)

def normalized_double_gaussian(x, mu1, sigma1, amp1, mu2, sigma2, amp2):
    """
    Computes a sum of two normalized Gaussian functions.
    """
    return (amp1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2) / np.sqrt(TWO_PI * sigma1 ** 2) +
            amp2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2) / np.sqrt(TWO_PI * sigma2 ** 2))

def continuum(x, a, b, c):
    """
    Represents a polynomial background function (continuum).
    """
    return a * x**2 + b * x + c

def combined_model_gaussian(x, a, b, c, mu, sigma, amplitude):
    """
    Combines a polynomial continuum with a single Gaussian peak model.
    """
    return continuum(x, a, b, c) + normalized_gaussian(x, mu, sigma, amplitude)

def combined_model_double_gaussian(x, a, b, c, mu1, sigma1, amp1, mu2, sigma2, amp2):
    """
    Combines a polynomial continuum with two Gaussian peaks.
    """
    return continuum(x, a, b, c) + normalized_double_gaussian(x, mu1, sigma1, amp1, mu2, sigma2, amp2)

def resolution_function(E, a, b, c):
    """
    Calculates the energy resolution as a function of energy.
    """
    return np.sqrt(a * E**-2 + b * E**-1 + c)

def linear_model(x, slope, intercept):
    """
    Linear model.
    """
    return slope * x + intercept