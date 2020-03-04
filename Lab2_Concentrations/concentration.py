import numpy as np
import math
import astropy
from scipy import stats
from scipy.integrate import quad
from astropy.stats import median_absolute_deviation

# TODO: Update comments

# Define some constants
point_a = -1 * math.sqrt(3)
point_b = math.sqrt(3)
scale_laplace = math.sqrt(2) / 2
scale_mixed = 3
sample_size = 100
mixed_size_a = 90
mixed_size_b = 10
trim = 0.05
upper_quartile_index = 74 # array started from index 0
lower_quartile_index = 24

# Samples: 
# 0 - Normal, 1 - Random, 2 - Laplace, 3 - Cauchy, 4 - Mixed
sample = [0.0 * sample_size] * 5

mean_square_dev      = [0] * 5 
average_absolute_dev = [0] * 5 
average_range        = [0] * 5 
inter_quartile_range = [0] * 5 # IQR
median_absolute_dev  = [0] * 5 # MAD

mean_square_dev_square       = [0] * 5
average_absolute_dev_square  = [0] * 5
average_range_square         = [0] * 5
inter_quartile_range_square  = [0] * 5
median_absolute_dev_square   = [0] * 5

def mean_square(sample):
    average = np.mean(sample)
    result = 0
    for i in range(len(sample)):
        result += pow(sample[i] - average, 2)
    return math.sqrt(result / len(sample))

def average_absolute(sample):
    median = np.median(sample)
    result = 0
    for i in range(len(sample)):
        result += abs(sample[i] - median)
    return result / len(sample)

def av_range(sample):
    # as long as sample is a sorted list
    return sample[len(sample)-1] - sample[0]

def iqr(sample):
    # as long as sample is a sorted list
    return sample[upper_quartile_index] - sample[lower_quartile_index] 


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

        print("\tFirst moment: {:.3f}".format(first_moment))
        print("\tSecond moment: {:.3f}".format(second_moment))
        print("\tDeviation: {:.3f}".format(deviation))

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
    sample[4] = np.concatenate([np.random.standard_normal(size = mixed_size_a), np.random.normal(0, scale_mixed, mixed_size_b)])

    for i in range(5):
        sample[i] = np.sort(sample[i])

    # Calculate all sums of characteristics of samples
    for i in range(5):
        mean_square_dev[i]      += mean_square(sample[i])
        average_absolute_dev[i] += average_absolute(sample[i])
        average_range[i]        += av_range(sample[i])
        inter_quartile_range[i] += stats.iqr(sample[i])
        median_absolute_dev[i]  += median_absolute_deviation(sample[i])

        mean_square_dev_square[i]      += pow(mean_square(sample[i]), 2)
        average_absolute_dev_square[i] += pow(average_absolute(sample[i]), 2)
        average_range_square[i]        += pow(av_range(sample[i]), 2)
        inter_quartile_range_square[i] += pow(stats.iqr(sample[i]), 2)
        median_absolute_dev_square[i]  += pow(median_absolute_deviation(sample[i]), 2)

print_mean_result('s', mean_square_dev, mean_square_dev_square)
print_mean_result('d', average_absolute_dev, average_absolute_dev_square)
print_mean_result('R', average_range, average_range_square)
print_mean_result('IQR', inter_quartile_range, inter_quartile_range_square)
print_mean_result('MAD', median_absolute_dev, median_absolute_dev_square)

#plt.hist(sample_mixed, bins=20)
#plt.savefig("mixed.pdf")


