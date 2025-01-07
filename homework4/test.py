import matplotlib.pyplot as plt
import numpy as np

file_list = [
    "./uniform_samples_100.dat",
    "./uniform_samples_1000.dat",
    "./uniform_samples_10000.dat",
    "./uniform_samples_100000.dat",
    "./normal_samples_100.dat",
    "./normal_samples_1000.dat",
    "./normal_samples_10000.dat",
    "./normal_samples_100000.dat",
]

for file in file_list:
    data = np.loadtxt(file)
    plt.hist(data, bins=100, alpha=0.5, label=file.split("/")[-1])
    plt.legend()
    plt.show()
