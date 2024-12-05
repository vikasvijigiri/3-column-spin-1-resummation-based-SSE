import os
import shutil
import erroranalysis as ea
import pandas as pd
import re
import numpy as np



# **********************************************************************************************************
#              GENERIC PART (NEVER CHANGE (except what you need exactly need is not found) 
# **********************************************************************************************************


class update_data:
    def __init__(self, mdir):
        # System parameters
        self.mdir = mdir
        print("Hello, you are extracting *.txt files from ", mdir)
        


    def extract_dat_files(self, root_folder, endswith):
        # Create the output folder if it doesn't exist
        # Create the output folder if it doesn't exist
        #os.makedirs(output_folder, exist_ok=True)

        #os.chmod(output_folder, 0o777)

        found_files_list = []
        # Iterate through all subdirectories and extract *.dat files
        for foldername, subfolders, filenames in os.walk(root_folder):
            #print(foldername,subfolders,filenames)
            for filename in filenames:
                if filename.endswith(endswith):
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


    # ******************************************************************

    def process_data(self, quoted_paths, process_extra_vars):
        data_list = []

        if process_extra_vars:
            # Instantiate extra_vars class
            ev = extra_vars()
        
        for file_path in quoted_paths:
            df = pd.read_csv(file_path, delim_whitespace=True)  
            num_lines = df.shape[0]   
            if num_lines == 0: continue         
            #nums = self.get_numbers_from_file(file_path)[:-depth]            
            if num_lines < 2:
                # Clone the row and append it as the next row
                dg = df.iloc[0].copy()
                df.loc[df.index.max() + 1] = dg
                #df = df.append(dg, ignore_index=True)
                print("\x1b[31m\"Warning:\"\x1b[0m"+ " Not enough data points for " + "\033[94m"+file_path+"\033[0m", " = ", num_lines)
                #print(df) 
            num_columns = df.shape[1]
            jackknife_results = [ea.jackknife(df.iloc[:, col].values) for col in range(num_columns)]

            if process_extra_vars:
                S1 = ev.param_arr(df, ["SMag_square_L"], [1])
                S2 = ev.param_arr(df, ["SMag_square_L"], [2])

                F1 = ev.param_arr(df, ["SMag_four_L"], [1])
                F2 = ev.param_arr(df, ["SMag_four_L"], [2])

                N1_corr = ev.param_arr(df, ["N_mx1_", "N_mx2_", "N_my1_", "N_my2_"], [0])
                V1_corr = ev.param_arr(df, ["V_mx1_", "V_mx2_", "V_my1_", "V_my2_"], [0])
                B1_corr = ev.param_arr(df, ["B_mx1_", "B_mx2_", "B_my1_", "B_my2_"], [0])

                N2_corr = ev.param_arr(df, ["N_mx1_", "N_mx2_", "N_my1_", "N_my2_"], [1])
                V2_corr = ev.param_arr(df, ["V_mx1_", "V_mx2_", "V_my1_", "V_my2_"], [1])
                B2_corr = ev.param_arr(df, ["B_mx1_", "B_mx2_", "B_my1_", "B_my2_"], [1])

                jackknife_results += [ea.jackknife_on_function(ev.SMag_binder, *np.concatenate((S1, F1)))]
                jackknife_results += [ea.jackknife_on_function(ev.SMag_binder, *np.concatenate((S2, F2)))]

                jackknife_results += [ea.jackknife_on_function(ev.Neel_order_params, *N1_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.Neel_order_params_ratio, *N1_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.VBS_order_params, *V1_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.VBS_order_params_ratio, *V1_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.VBS_order_params, *B1_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.VBS_order_params_ratio, *B1_corr)]

                jackknife_results += [ea.jackknife_on_function(ev.Neel_order_params, *N2_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.Neel_order_params_ratio, *N2_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.VBS_order_params, *V2_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.VBS_order_params_ratio, *V2_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.VBS_order_params, *B2_corr)]
                jackknife_results += [ea.jackknife_on_function(ev.VBS_order_params_ratio, *B2_corr)]

            # jackknife_results += [ea.jackknife_on_function(self.Smag_binder, df.iloc[:, 1].values, df.iloc[:, 2].values)]
            #binder_result = ea.jackknife_on_function(Smag_binder, df.iloc[:, 1].values, df.iloc[:, 2].values)
            darr = [item for sublist in jackknife_results for item in sublist] #+ list(binder_result)
            data_list.append(darr)  
        return data_list
        #    return df.values  # Assuming you want to return the entire DataFrame if there's not enough data





    def avg_over_cols_inside_file(self, found_files, process_extra_vars=False):
        quoted_paths = [f"{os.path.normpath(path)}" for path in found_files]
        #print(quoted_paths)
        f = quoted_paths[0]
        df = pd.read_csv(f, delim_whitespace=True) 
        headers = df.columns.tolist()
        if process_extra_vars: 
            headers += ["SMag_binder_L1", "SMag_binder_L2", "O_N_L1", "R_N_L1", "O_V_L1", "R_V_L1", 
                        "O_B_L1", "R_B_L1", "O_N_L2", "R_N_L2", "O_V_L2", "R_V_L2", 
                        "O_B_L2", "R_B_L2"]
        headers = [f'{item}{suffix}' for item in headers for suffix in ["_avg", "_err"]]
        #print(headers)
        data_list = self.process_data(quoted_paths, process_extra_vars)
        return pd.DataFrame(data_list, columns=headers)
                                 


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
                else:
                    print(f"Variable '{pattern}' is not found in naming the file. Either remove it or use naming convention that involves this parameter")
            #print(data)
            df[rf'{variables[i]}'] =  data 
            df = df.apply(pd.to_numeric, errors='coerce')
 
        return df        


# **********************************************************************************************************
#                 Model specific processing part (YOU CAN ADD YOUR OWN HERE) 
# **********************************************************************************************************

class extra_vars():
    def __init__(self):
        # Global parameters (if any)
        print("You are now processing observables outside RAW data files.")


    # ******************** Observables part **************************

#    def SMag_square(self, df, sym=True):
#        if sym: SMag_square_kernel = np.array([1.,1.,1.])
#        else: SMag_square_kernel = np.array([1.,1.,-1.])
#        N = []
#        N.append(df["SMag_squareL1"]*SMag_square_kernel[0])
#        N.append(df["SMag_squareL2"]*SMag_square_kernel[1]) 
#        N.append(df["SMag_squareI1"]*SMag_square_kernel[2])
#        return N


#    def SMag_four(self, df, sym=True):
#        if sym: SMag_square_kernel = np.array([1.,1.,1.,1.,1.])
#        else: SMag_square_kernel = np.array([1.,1.,1.,-1,-1])
#        N = []
#        N.append(df["SMag_fourL1"]*SMag_square_kernel[0])
#        N.append(df["SMag_fourL2"]*SMag_square_kernel[1]) 
#        N.append(df["SMag_fourI1"]*SMag_square_kernel[2])
#        N.append(df["SMag_fourI2"]*SMag_square_kernel[3])
#        N.append(df["SMag_fourI3"]*SMag_square_kernel[4]) 
#        return N

    def param_arr(self, df, param_type, layers):
        # Initialize an empty list to store the result variables
        result_variables = []
        # param_type = ["N_mx_", "N_my_"]
        # Loop through the indices
        for par in param_type:
            for i in layers:
                #variable_name = f"N{i * 4 + j + 1}"
                result_variable = df[f"{par}{i}"].values 
                #locals()[variable_name] = result_variable
                result_variables.append(result_variable)  
        return result_variables


    #*******************************************************************************************


    def SMag_binder(self, *SMag_square_four_arr):
        NS1, NF1 = SMag_square_four_arr
        return  2.5*(1. - (NF1)/ (3. * (NS1)**2))

    def VBS_order_params_ratio(self, *Nreal_arr):
        N1, N2, N3, N4 = Nreal_arr
        return 1. - (N1 / N2)


    def VBS_order_params(self, *Nreal_arr):
        N1, N2, N3, N4 = Nreal_arr   
        return N2

    def Neel_order_params_ratio(self, *Nreal_arr):
        N1, N2, N3, N4 = Nreal_arr
        return 1. - (N1 / N2 + N3 / N4) / 2.


    def Neel_order_params(self, *Nreal_arr):
        N1, N2, N3, N4 = Nreal_arr   
        return (N2 + N4) / 2.

    # *****************************************************************************************



# *********************************************************************************************
#                             MAIN (Access)
# *********************************************************************************************



def run_update(root_folder, variables, endswith='data_avg.txt', process_extra_obs=False):  
    upd = update_data(root_folder)      
    found_files = upd.extract_dat_files(root_folder, endswith)


#    if process_extra_obs: 
#        pro = extra_vars()

    try:
        # Check if the file is empty
        if len(found_files) < 1:
            raise ValueError(f"The file list is empty.")
        
        df  = upd.avg_over_cols_inside_file(found_files, process_extra_obs)
        df_match    = upd.match_pattern_file(found_files, variables)
        df = pd.concat([df_match, df], axis=1)
        df = df.sort_values(by=variables[:], ascending=True)
        df = df.reset_index(drop=True)
        
        #if endswith != 'data_avg.txt': df = df.transpose().reset_index(drop=True)

        # Return some result if needed
        return found_files, df

    except ValueError as e:
        # Handle the case where the file is empty
        print(f"Error: {e}")
        # Return some error message or value
        return [], pd.DataFrame()



def savefile(dg, output_dir, ofile_path):
    os.makedirs(output_dir, exist_ok=True)
    os.chmod(output_dir, 0o777)
    file_path = output_dir + ofile_path
    dg.to_string(file_path, sparsify=False, index=None, float_format='%.8e')
    

