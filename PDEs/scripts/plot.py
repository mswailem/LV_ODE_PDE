import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from multiprocessing import Pool, Value, Lock

# This function finds the global minimum and maximum values for u and v in even and odd frames seperately
def find_global_min_max(files):
    even_min_u, even_max_u = float('inf'), float('-inf')
    even_min_v, even_max_v = float('inf'), float('-inf')
    odd_min_u, odd_max_u = float('inf'), float('-inf')
    odd_min_v, odd_max_v = float('inf'), float('-inf')
    
    for file in files:
        # Extract time step from the file name
        t = int(file.replace('_x.dat', ''))
        
        df = pd.read_csv(file, sep=r'\s+', header=None)
        u, v = df[2], df[3]
        
        if t % 2 == 0:  # Even frame
            even_min_u = min(even_min_u, u.min())
            even_max_u = max(even_max_u, u.max())
            even_min_v = min(even_min_v, v.min())
            even_max_v = max(even_max_v, v.max())
        else:  # Odd frame
            odd_min_u = min(odd_min_u, u.min())
            odd_max_u = max(odd_max_u, u.max())
            odd_min_v = min(odd_min_v, v.min())
            odd_max_v = max(odd_max_v, v.max())
    
    return even_min_u, even_max_u, even_min_v, even_max_v, odd_min_u, odd_max_u, odd_min_v, odd_max_v

# Function that processes a single file (a single timestep)
def process_file(file):
    # TODO: Get rid of global variables, keeping this for now as it works, maybe this was needed for parallel processing?
    global progress_counter, lock, num_files
    global even_min_u, even_max_u, even_min_v, even_max_v
    global odd_min_u, odd_max_u, odd_min_v, odd_max_v

    # Extract time step from the file name
    t = int(file.replace('_x.dat', ''))

    # Load data from the file
    df = pd.read_csv(file, sep=r'\s+', header=None)
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

    # Interpolate red and blue intensities just so that the image is not pixelated
    Z_red = griddata(points, red_intensity, (X, Y), method='cubic')
    Z_blue = griddata(points, blue_intensity, (X, Y), method='cubic')

    # Normalize using global min and max values for even and odd frames (NOTE: This was done for a specific dataset, and may not work as nicely for other datasets)
    if t % 2 == 0:  # Even frame
        Z_blue = (Z_blue - even_min_v) / (even_max_v - even_min_v)
    else:  # Odd frame
        Z_red = (Z_red - odd_min_u) / (odd_max_u - odd_min_u)

    # Combine red and blue intensities into an RGB array
    RGB = np.dstack((Z_red, np.zeros_like(Z_red), Z_blue))

    # Plot as an image
    plt.figure(figsize=(8, 8))
    plt.imshow(RGB, origin='lower', extent=(min(x), max(x), min(y), max(y)))

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
    file_path = os.path.join('../images', file)
    if os.path.isfile(file_path):
        os.remove(file_path)

# Get all dataset files
dataset_files = [file for file in sorted(os.listdir('.')) if file.endswith('_x.dat')]

# Update global min and max values
even_min_u, even_max_u, even_min_v, even_max_v, odd_min_u, odd_max_u, odd_min_v, odd_max_v = find_global_min_max(dataset_files)

# Global variables for tracking progress
progress_counter = Value('i', 0)
num_files = len(dataset_files)
lock = Lock()

# Parallel processing using multiprocessing.Pool
with Pool() as pool:
    pool.map(process_file, dataset_files)

print("\nPlotting complete.")
