# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from natsort import natsorted
# %%
path = "."
output_path = "../processed_output"
os.makedirs(output_path, exist_ok=True)

file_list = natsorted(glob.glob(os.path.join(path, "*_x.dat")))

std_values = []

for file in file_list:
    data = pd.read_csv(file, sep=r'\s+', header=None)
    data.columns = ["x", "y", "u", "v"]
    std = data["u"].std()
    time = int(os.path.basename(file).split('_')[0])
    std_values.append({"time": time, "std": std})

result = pd.DataFrame(std_values)
result.to_csv(os.path.join(output_path, "std.dat"), index=False)
