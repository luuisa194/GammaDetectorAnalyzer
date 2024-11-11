
from numpy import inf as INF
import numpy as np

#the following functions are from the provided documents

def first_moment(x, y):
    return np.sum(x * y) / np.sum(y)

def second_moment(x, y):
    x0 = first_moment(x, y)
    return np.sum((x - x0) ** 2 * y) / np.sum(y)

def gaussian_initial_estimates(channels, counts):
    """Estimates the three parameters of a Gaussian distribution."""
    channels = np.asarray(channels)  # Convert channels to a numpy array
    counts = np.asarray(counts)      # Convert counts to a numpy array
    
    mu0 = first_moment(channels, counts)
    sig0 = np.sqrt(second_moment(channels, counts))
    a0 = np.sum(counts)
    return mu0, sig0, a0

def in_interval(x, xmin=-INF, xmax=INF):
    """Boolean mask with value True for x in [xmin, xmax)."""
    _x = np.asarray(x)
    return np.logical_and(xmin <= _x, _x < xmax)

def filter_in_interval(x, y, xmin, xmax):
    """Selects only elements of x and y where xmin <= x < xmax."""
    # Ensure x and y have the same shape
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")
    
    # Create mask and apply to both arrays
    _mask = in_interval(x, xmin, xmax)
    return [np.asarray(arr)[_mask] for arr in (x, y)]
