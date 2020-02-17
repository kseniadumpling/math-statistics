import numpy as np
import math
from scipy import stats
#from matplotlib import pyplot as plt
#import mixture

# Define some constants
point_a = -1 * math.sqrt(3)
point_b = math.sqrt(3)
scale_laplace = math.sqrt(2) / 2
scale_mixed = 3
sample_size = 100
mixed_size_a = 90
mixed_size_b = 10
trim = 0.05

# Samples: 
# 0 - Normal, 1 - Random, 2 - Laplace, 3 - Cauchy, 4 - Mixed
sample = [0.0 * sample_size] * 5

mean    = [0] * 5
median  = [0] * 5
extreme = [0] * 5
trimmed = [0] * 5

mean_square    = [0] * 5
median_square  = [0] * 5
extreme_square = [0] * 5
trimmed_square = [0] * 5

# Finding half sum of extreme elements
def extr_mean(sample):
    sample = np.sort(sample)
    min = sample[0]
    max = sample[sample_size-1]
    return (min + max) / 2 

# Printing results of calculations
# Input: 
#   name (string) - name of characteristic (mean
#   mean (float) - sum of 1000 calculations
#   mean_square (float) - sum of squares of 1000 calculations
def print_mean_result(name, mean, mean_square):
    print("\n\n ---------- {:s} ----------\n".format(name))
    for i in range(5):
        print("\n{:d} distribution \n".format(i))

        first_moment = mean[i] / 1000
        second_moment = mean_square[i] / 1000
        deviation = second_moment - pow(first_moment, 2)

        print("\tFirst moment: {:.10f}".format(first_moment))
        print("\tSecond moment: {:.10f}".format(second_moment))
        print("\tDeviation: {:.10f}".format(deviation))

# ----- Main loop -----
# Monte Carlo method
# The mean, median, extremal mean and trimmed mean will be calculated 
# one thousand times and after there will be calculated average of first & second moments
for iteration in range(1000):
    # Generate samples of different kinds of distributions
    sample[0] = np.random.standard_normal(size = sample_size)
    sample[1] = (point_b - point_a) * np.random.random_sample(size = sample_size) + point_a
    sample[2] = np.random.laplace(0, scale_laplace, sample_size)
    sample[3] = np.random.standard_cauchy(size = sample_size)
    sample[4]  = np.concatenate([np.random.standard_normal(size = mixed_size_a), np.random.normal(0, scale_mixed, mixed_size_b)])

    # Calculate all sums of characteristics of samples
    for i in range(5):
        mean[i]    += np.mean(sample[i])
        median[i]  += np.median(sample[i])
        extreme[i] += extr_mean(sample[i])
        trimmed[i] += stats.trim_mean(sample[i], trim)

        mean_square[i]    += pow(np.mean(sample[i]), 2)
        median_square[i]  += pow(np.median(sample[i]), 2)
        extreme_square[i] += pow(extr_mean(sample[i]), 2)
        trimmed_square[i] += pow(stats.trim_mean(sample[i], trim), 2)

print_mean_result('Mean', mean, mean_square)
print_mean_result('Median', median, median_square)
print_mean_result('Extreme', extreme, extreme_square)
print_mean_result('Trimmed', trimmed, trimmed_square)

#plt.hist(sample_mixed, bins=20)
#plt.savefig("mixed.pdf")


