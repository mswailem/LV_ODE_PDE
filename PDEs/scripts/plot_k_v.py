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
    t = file.replace('_k.dat', '')

    # Load data from the file
    df = pd.read_csv(file, sep=r'\s+', header=None)
    x, y = df[0], df[1]
    v = df[3]

    # Prepare coordinates and intensity
    points = np.column_stack((x, y))
    intensity = v

    # Create a grid for interpolation
    xi = np.linspace(min(x), max(x), 128)
    yi = np.linspace(min(y), max(y), 128)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate intensity
    Z = griddata(points, intensity, (X, Y), method='nearest')

    # Remove DC signal at wave vector (0,0)
    center_x = np.argmin(np.abs(xi))
    center_y = np.argmin(np.abs(yi))
    Z[center_y, center_x] = 0

    # Normalize intensity using a logarithmic scale to improve contrast
    Z_log = np.log1p(Z)  # log(1 + Z) to avoid issues with log(0)
    Z_normalized = Z_log / Z_log.max()

    # Create image using a perceptually uniform colormap
    plt.figure(figsize=(8, 8))
    im = plt.imshow(Z_normalized, origin='lower', extent=(-0.05, 0.05, -0.05, 0.05), cmap='gray')

    # Add a colorbar
    cbar = plt.colorbar(im, ax=plt.gca())
    cbar.set_label('Normalized Intensity (log scale)')

    # Add a title
    plt.title(f't={t}')

    # Save the plot
    plt.savefig(f'../images_k_v/{t}.png')
    plt.close()

    # Update progress counter
    with lock:
        progress_counter.value += 1
        print(f"Processing time step: {progress_counter.value}/{num_files}", end='\r', flush=True)

# Set up directories
os.makedirs('../images_k_v', exist_ok=True)
for file in os.listdir('../images_k_v'):
    file_path = os.path.join('../images_k_v', file)
    if os.path.isfile(file_path):  # Check if the item is a file
        os.remove(file_path)

# Get all dataset files
dataset_files = [file for file in sorted(os.listdir('.')) if file.endswith('_k.dat')]

# Global variables for tracking progress
progress_counter = Value('i', 0)
num_files = len(dataset_files)
lock = Lock()

# Parallel processing using multiprocessing.Pool
with Pool() as pool:
    pool.map(process_file, dataset_files)

print("\nPlotting complete.")
