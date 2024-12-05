import numpy as np
import os
import shutil
import erroranalysis as ea
import pandas as pd
import re





class get_files:
    def __init__(self, mdir):
        # System parameters
        self.mdir = mdir
        #print("Hello, you are extracting *series.txt files from ", mdir)
        



    def extract_dat_files(self, root_folder):
        # Create the output folder if it doesn't exist
        # Create the output folder if it doesn't exist
        #os.makedirs(output_folder, exist_ok=True)

        #os.chmod(output_folder, 0o777)

        found_files_list = []
        # Iterate through all subdirectories and extract *.dat files
        for foldername, subfolders, filenames in os.walk(root_folder):
            #print(foldername,subfolders,filenames)
            for filename in filenames:
                if filename.endswith('series.txt'):
                    file_path = os.path.join(foldername, filename)
                    #print(f"Found: {file_path}")
                    depth = foldername.count(os.path.sep) - root_folder.count(os.path.sep)

                    found_files_list.append(file_path)                    

                    # Copy the *.dat file to the output folder
                    #shutil.copy(file_path, output_folder)
                    #print(f"Extracted: {file_path} to {output_folder}/{filename}")
        return found_files_list

    def print_found_files(self, found_files_list): 
        # Now found_files_list contains all the file paths with the specified extension
        print("List of found files:")
        for file_path in found_files_list:
            print(file_path)

    def match_pattern_file(self, file_path_list, variables):
        l = len(variables)
        patterns = [rf'{variables[i]}.*?(\d+(\.\d+)?)' for i in range(l)]

        df = pd.DataFrame()
        # Iterate over patterns
        for i, pattern in enumerate(patterns):  
            data = []
            for file_path in file_path_list:
                match = re.search(pattern, file_path)

                if match:
                    number_right_of_pattern = match.group(1)
                    data.append(number_right_of_pattern)
                    #print(f"Number right of the pattern '{pattern}': {number_right_of_pattern}")
                #else:
                #    print(f"Variable '{pattern}' is not found in naming the file. Either remove it or use naming convention that involves this parameter")
            #print(data)
            df[rf'{variables[i]}'] =  data 
            df = df.apply(pd.to_numeric, errors='coerce')
 
        return df   

class autocorrelation:
    def __init__(self, mdir):
        # System parameters
        self.mdir = mdir
        #print("Hello, you are analyzing time_series files from ", mdir)
        

    def auto_correlation(self, data, lag):
        n = len(data)
        mean = np.mean(data)
        numerator = np.sum((data[:n-lag] - mean) * (data[lag:] - mean))
        denominator = np.sum((data - mean)**2)
        return numerator / denominator if denominator != 0 else 0

    def calculate_integrated_auto_correlation_time(self, data, max_lag):
        integrated_auto_correlation_time = 0.0  # Initial value
        #print(data)
        for lag in range(0, max_lag + 1):
            correlation = self.auto_correlation(data, lag)
            if correlation > 0:
                integrated_auto_correlation_time += correlation
            else:
                #print("Error: Tau is negative!")
                break  # Stop if correlation becomes negative

        return integrated_auto_correlation_time

    def calculate_auto_correlation_time(data, max_lag):
        auto_correlation_times = []
        for lag in range(1, max_lag + 1):
            correlation = self.auto_correlation(data, lag)
            if correlation > 0:
                auto_correlation_times.append(lag * correlation)
            else:
                break  # Stop if correlation becomes negative

        return np.sum(auto_correlation_times)

    def block_resampling(self,data, block_size):
        n = len(data)
        num_blocks = n // block_size
        block_means = np.mean(data[:num_blocks * block_size].reshape((num_blocks, block_size)), axis=1)
        return block_means

    def estimate_error(self,data, max_lag, num_blocks):
        n = len(data)
        block_sizes = [2**k for k in range(int(np.log2(n)) - 1)]  # Using 5 different block sizes
        auto_corr_times = []

        for block_size in block_sizes:
            block_means = self.block_resampling(data, block_size)
            auto_corr_time = self.calculate_integrated_auto_correlation_time(block_means, max_lag)
            auto_corr_times.append(auto_corr_time)

        mean_auto_corr_time = np.mean(auto_corr_times)
        std_dev_auto_corr_time = np.std(auto_corr_times)
        
        # Standard error of the mean
        error = std_dev_auto_corr_time / np.sqrt(num_blocks)
        return mean_auto_corr_time, error


def run_autocorr_times(root_folder, variables, obs_compute_tau):  
    upd = get_files(root_folder)      
    found_files = upd.extract_dat_files(root_folder)
    df_match = upd.match_pattern_file(found_files, variables)
    #upd.print_found_files(found_files) 
   
    ac = autocorrelation(root_folder)
#    for obs in obs_compute_tau:
    obs_tau = [f'{item}{suffix}' for item in obs_compute_tau for suffix in ["_tau_avg", "_tau_err"]]
    df_tau   = pd.DataFrame([], columns=obs_tau) 
    df = pd.concat([df_match, df_tau], axis=1)   


    for i, obs in enumerate(obs_compute_tau):
        tau_mean = []; tau_err = []        
        for ifile in found_files:    
            data = pd.read_csv(ifile, delim_whitespace=True, comment='#', usecols = [obs])
            data = np.asarray(data).ravel()    
            auto_corr = ac.calculate_integrated_auto_correlation_time(data, 15)
            T_mean, T_err = ac.estimate_error(data, 20, 10000)
            tau_mean.append(T_mean)
            tau_err.append(T_err)
        df[obs_tau[2*i]] = tau_mean
        df[obs_tau[2*i+1]] = tau_err
#    df_tau   = pd.DataFrame(tau_list, columns=['Tau'])
#    df = pd.concat([df_match, df_tau], axis=1)
    df = df.sort_values(by='L', ascending=True)
    df = df.reset_index(drop=True)
    return df 


