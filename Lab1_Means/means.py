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

# ----- Main loop -----
# Monte Carlo method
# The mean, median, extremal mean and trimmed mean will be computed 
# one thousand times and after there will be average first & second moments
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

# Compute mean moments 
print("\n\n ---------- Mean ----------\n")
for i in range(5):
    print("\n{:d} distribution \n".format(i))

    first_moment = mean[i] / 1000
    second_moment = mean_square[i] / 1000
    deviation = second_moment - pow(first_moment, 2)

    print("\tFirst moment: {:.10f}".format(first_moment))
    print("\tSecond moment: {:.10f}".format(second_moment))
    print("\tDeviation: {:.10f}".format(deviation))

# Compute median moments 
print("\n\n ---------- Median ----------\n")
for i in range(5):
    print("\n{:d} distribution \n".format(i))

    first_moment = median[i] / 1000
    second_moment = median_square[i] / 1000
    deviation = second_moment - pow(first_moment, 2)

    print("\tFirst moment: {:.10f}".format(first_moment))
    print("\tSecond moment: {:.10f}".format(second_moment))
    print("\tDeviation: {:.10f}".format(deviation))


# Compute half sum of extreme elements moments 
print("\n\n ---------- Extreme ----------\n")
for i in range(5):
    print("\n{:d} distribution \n".format(i))

    first_moment = extreme[i] / 1000
    second_moment = extreme_square[i] / 1000
    deviation = second_moment - pow(first_moment, 2)

    print("\tFirst moment: {:.10f}".format(first_moment))
    print("\tSecond moment: {:.10f}".format(second_moment))
    print("\tDeviation: {:.10f}".format(deviation))

# Compute trimmed mean moments 
print("\n\n ---------- Trimmed ----------\n")
for i in range(5):
    print("\n{:d} distribution \n".format(i))

    first_moment = trimmed[i] / 1000
    second_moment = trimmed_square[i] / 1000
    deviation = second_moment - pow(first_moment, 2)

    print("\tFirst moment: {:.10f}".format(first_moment))
    print("\tSecond moment: {:.10f}".format(second_moment))
    print("\tDeviation: {:.10f}".format(deviation))


#plt.hist(sample_mixed, bins=20)
#plt.savefig("mixed.pdf")


