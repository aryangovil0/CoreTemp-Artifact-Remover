import os
import numpy as np
import pandas as pd
from scipy import interpolate
import re

def process_file(file_path):
    df = pd.read_excel(file_path)
    
    # Find the first row that can be parsed as a datetime
    start_row = 0
    for i, value in enumerate(df.iloc[:, 0]):
        try:
            pd.to_datetime(value)
            start_row = i
            break
        except ValueError:
            continue
    
    # If no valid datetime found, raise an error
    if start_row == len(df):
        raise ValueError(f"No valid datetime found in the first column of {file_path}")
    
    # Use only the rows from start_row onwards
    df = df.iloc[start_row:].reset_index(drop=True)
    
    temp0 = df.iloc[:, -1].values
    datetime_col = pd.to_datetime(df.iloc[:, 0])
    HH = datetime_col.dt.hour.to_numpy()
    MM = datetime_col.dt.minute.to_numpy()

    # Handle hour jumps due to midnight resets
    JUMP = 0
    for kk in range(1, len(HH)):
        if HH[kk] == 0 and MM[kk] == 0:
            JUMP += 24
        HH[kk] += JUMP

    # Generate uniform time series
    time0 = (HH - HH[0]) * 60 + MM
    time = np.arange(time0[-1] - time0[0] + 1) + time0[0]

    # Map raw temperatures to uniform time series
    temp1 = np.full(len(time), np.nan)
    temp1[time0 - time0[0]] = temp0
    idx_nan = np.isnan(temp1)
    valid_idx = ~np.isnan(temp1)

    # Interpolate missing values
    f = interpolate.interp1d(
        np.where(valid_idx)[0], 
        temp1[valid_idx], 
        kind='linear', 
        bounds_error=False, 
        fill_value="extrapolate"
    )
    temp1 = f(np.arange(len(temp1)))

    # Apply invalidation rules
    idx_out_of_range = np.where((temp1 <= 35.5) | (temp1 >= 39))[0]
    dtemp = np.diff(temp1) / np.diff(time)
    idx_abrupt_change = np.union1d(
        np.where(dtemp < -0.4/3)[0], 
        np.where(dtemp > 0.4/3)[0]
    )
    idx_invalid = np.union1d(idx_out_of_range, idx_abrupt_change)
    temp1[idx_invalid] = np.nan

    # Identify and handle invalid segments
    invalid_segments = ''.join(map(str, np.isnan(temp1).astype(int)))
    segment_lengths = [len(m.group()) for m in re.finditer('1+', invalid_segments)]
    segment_starts = [m.start() for m in re.finditer('1+', invalid_segments)]
    for kk, seg_len in enumerate(segment_lengths):
        if seg_len >= 1:
            idx_invalid = np.union1d(
                idx_invalid, 
                np.arange(max(0, segment_starts[kk]-10), min(len(temp1), segment_starts[kk]+seg_len+20))
            )

    temp = temp1.copy()
    temp[idx_invalid] = np.nan

    # Determine start time
    if np.isnan(temp[0]):
        first_valid = np.where(~np.isnan(temp))[0][0]
        Starttime = datetime_col.iloc[first_valid].time()
        temp = temp[first_valid:]
        time = time[first_valid:]
    else:
        Starttime = datetime_col.iloc[0].time()

    return temp, time, Starttime, temp0, time0

def impute(x0, bdry, LIB, LIBmean):
    x = x0.copy()
    idx_nan = np.isnan(x)

    if not idx_nan.any():
        return x

    # Locate missing segments
    diff_idx_nan = np.diff(np.where(idx_nan)[0])
    break_points = np.where(diff_idx_nan > 1)[0]
    seg_start = np.concatenate(([np.where(idx_nan)[0][0]], np.where(idx_nan)[0][break_points + 1]))
    seg_end = np.concatenate((np.where(idx_nan)[0][break_points], [np.where(idx_nan)[0][-1]]))

    for start, end in zip(seg_start, seg_end):
        if start == 0:
            neighbors = x[end+1:min(end+1+bdry, len(x))]
        elif end == len(x) - 1:
            neighbors = x[max(0, start-bdry):start]
        else:
            neighbors = np.concatenate((x[max(0, start-bdry):start], x[end+1:min(end+1+bdry, len(x))]))

        if len(neighbors) == 0:
            # If no neighbors are found, skip imputation for this segment
            continue

        x_mean = np.mean(neighbors)
        best_match, min_error = None, np.inf
        for lib, lib_mean in zip(LIB, LIBmean):
            if len(lib) >= end - start + 1:
                lib_segment = lib[:end-start+1] + lib_mean
                error = np.linalg.norm(lib_segment - x_mean)
                if error < min_error:
                    best_match, min_error = lib_segment, error

        if best_match is None:
            # If no match is found in the library, use linear interpolation
            idx_valid = ~np.isnan(x)
            f = interpolate.interp1d(
                np.where(idx_valid)[0], 
                x[idx_valid], 
                kind='linear', 
                bounds_error=False, 
                fill_value="extrapolate"
            )
            x[start:end+1] = f(np.arange(start, end+1))
        else:
            # Assign the best match
            x[start:end+1] = best_match

    # Interpolate remaining NaNs if any
    idx_valid = ~np.isnan(x)
    f = interpolate.interp1d(
        np.where(idx_valid)[0], 
        x[idx_valid], 
        kind='linear', 
        bounds_error=False, 
        fill_value="extrapolate"
    )
    return f(np.arange(len(x)))

def main(folder_path):
    files = [f for f in os.listdir(folder_path) if f.endswith('.xlsx')]
    output_folder = os.path.join(folder_path, 'output_data_test')
    os.makedirs(output_folder, exist_ok=True)

    LIB, LIBmean = [], []

    for file in files:
        print(f"Processing {file}...")
        file_path = os.path.join(folder_path, file)
        temp, time, Starttime, temp0, time0 = process_file(file_path)

        # Update library with valid segments
        valid_segments = ''.join(map(str, ~np.isnan(temp).astype(int)))
        segment_lengths = [len(m.group()) for m in re.finditer('1+', valid_segments)]
        segment_starts = [m.start() for m in re.finditer('1+', valid_segments)]
        for start, length in zip(segment_starts, segment_lengths):
            segment = temp[start:start+length] - np.mean(temp[~np.isnan(temp)])
            LIB.append(segment)
            LIBmean.append(np.mean(temp[~np.isnan(temp)]))

        # Perform imputation
        imputed_temp = impute(temp, bdry=60, LIB=LIB, LIBmean=LIBmean)

        # Save the output
        save_output(output_folder, file, temp, time, Starttime, temp0, time0, imputed_temp)

    print(f"All files processed. Results saved in {output_folder}")

def save_output(output_folder, file, temp, time, Starttime, temp0, time0, imputed_temp):
    base_name = os.path.splitext(file)[0]
    
    # Save raw and processed data
    raw_data_path = os.path.join(output_folder, f"{base_name}_raw.csv")
    processed_data_path = os.path.join(output_folder, f"{base_name}_processed.csv")
    imputed_data_path = os.path.join(output_folder, f"{base_name}_imputed.csv")
    
    # Save Starttime
    starttime_path = os.path.join(output_folder, f"{base_name}_starttime.txt")
    with open(starttime_path, 'w') as f:
        f.write(str(Starttime))
    
    # Save raw data
    pd.DataFrame({
        "Time (raw)": time0,
        "Temperature (raw)": temp0
    }).to_csv(raw_data_path, index=False)
    
    # Save processed data
    pd.DataFrame({
        "Time (processed)": time,
        "Temperature (processed)": temp
    }).to_csv(processed_data_path, index=False)
    
    # Save imputed data
    pd.DataFrame({
        "Time (imputed)": time,
        "Temperature (imputed)": imputed_temp
    }).to_csv(imputed_data_path, index=False)

if __name__ == "__main__":
    folder_path = ""
    main(folder_path)
