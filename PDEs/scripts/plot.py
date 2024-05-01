import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from multiprocessing import Pool, Value, Lock

# Helper function to generate a plot for a single dataset file
def process_file(file):
    # Global variables for progress tracking
    global progress_counter, lock, num_files

    # Extract time step from the file name
    t = file.replace('_x.dat', '')

    # Load data from the file
    df = pd.read_csv(file, delim_whitespace=True, header=None)
    x, y = df[0], df[1]
    u, v = df[2], df[3]

    # Prepare coordinates and colors
    points = np.column_stack((x, y))
    red_intensity = u
    blue_intensity = v

    # Create a grid for interpolation
    xi = np.linspace(min(x), max(x), 500)
    yi = np.linspace(min(y), max(y), 500)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate red and blue intensities
    Z_red = griddata(points, red_intensity, (X, Y), method='cubic')
    Z_blue = griddata(points, blue_intensity, (X, Y), method='cubic')

    # Combine red and blue intensities into an RGB array
    RGB = np.dstack((Z_red / Z_red.max(), np.zeros_like(Z_red), Z_blue / Z_blue.max()))

    # Plot as an image
    plt.figure(figsize=(8, 8))
    plt.imshow(RGB, origin='lower', extent=[min(x), max(x), min(y), max(y)])

    # Add a title
    plt.title(f't={t}')

    # Save the plot
    plt.savefig(f'../images/{t}.png')
    plt.close()

    # Update progress counter
    with lock:
        progress_counter.value += 1
        print(f"Processing time step: {progress_counter.value}/{num_files}", end='\r', flush=True)

# Set up directories
os.makedirs('../images', exist_ok=True)
for file in os.listdir('../images'):
    os.remove(f'../images/{file}')

# Get all dataset files
dataset_files = [file for file in sorted(os.listdir('.')) if file.endswith('_x.dat')]

# Global variables for tracking progress
progress_counter = Value('i', 0)
num_files = len(dataset_files)
lock = Lock()

# Parallel processing using multiprocessing.Pool
with Pool() as pool:
    pool.map(process_file, dataset_files)

print("\nPlotting complete.")

