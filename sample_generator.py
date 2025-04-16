import numpy as np
import pandas as pd

# Function to generate synthetic data
def generate_synthetic_data(size=1000):
    # Simulate 'auc' as a random value between 0 and 1 (e.g., area under the curve)
    auc = np.random.uniform(0, 1, size)
    
    # Simulate 'length' as a random value between 100 and 1000 (e.g., genomic region length)
    length = np.random.lognormal(mean=3, sigma=1, size=size).astype(int)  # Log-normal for length distribution
    
    # Combine into a DataFrame
    data = pd.DataFrame({
        'auc': auc,
        'length': length
    })
    
    return data

# Function to save data with overlap between two files
def save_data_with_overlap(size):
    print(f"Started working on: Sample size: {size}")
    # Generate two sets of data
    data1 = generate_synthetic_data(size)
    data2 = generate_synthetic_data(size)
    
    # Define the overlap size (e.g., 50% overlap)
    overlap_size = size // 2
    
    # Take the first 'overlap_size' rows from data1 and merge them with the rest of data2
    overlap_data = data1.head(overlap_size)
    unique_data2 = data2.tail(size - overlap_size)
    
    # Merge the two datasets: overlap + unique part of data2
    data2_with_overlap = pd.concat([overlap_data, unique_data2]).reset_index(drop=True)
    
    # Save the files
    file_name_1 = f'synthetic_data_1_{size}.auc'
    file_name_2 = f'synthetic_data_2_{size}.auc'
    
    data1.to_csv(file_name_1,sep='\t', index=False, header=False, float_format="%.6f")
    data2_with_overlap.to_csv(file_name_2, sep='\t', index=False, header=False, float_format="%.6f")
    
    print(f"Generated files: {file_name_1} and {file_name_2}")

# Generate and save data for the desired sizes
sizes = [10000000]  # Sizes you want to simulate

for size in sizes:
    save_data_with_overlap(size)

