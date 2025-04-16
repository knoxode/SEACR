import argparse
import time
from pathlib import Path
import numpy as np
from scipy.stats import gaussian_kde
from statsmodels.distributions import ECDF as ecdf

#Parsing input arguments
parser = argparse.ArgumentParser(
        prog="pySEACR",
        description="Calculate area under the curve threshold for CUT&RUN peaks. Implemented in Python and NumPy"
        )

parser.add_argument('-E', '--exp', type= lambda p: Path(p).resolve(), help="Input AUC values from experiment CUT&RUN - must be a file", required=True)
parser.add_argument('-C', '--ctrl', help="Input AUC values from control CUT&RUN - must be a file or number", required=True)
parser.add_argument('-N',"--norm", choices=['yes', 'no'], help="Whether to normalize control and experimental files", required=True)
parser.add_argument('-O', "--output", help="Output Prefix")

args = parser.parse_args()

########################
#  CORE FUNCTIONALITY  #
########################

def dist2d(x, y, x0=0, y0=0):
    dx, dy = 1 - x0, 1 - y0
    mag = np.hypot(dx, dy)
    dx, dy = dx / mag, dy / mag
    px, py = x - x0, y - y0
    proj = px * dx + py * dy
    closest_x = x0 + proj * dx
    closest_y = y0 + proj * dy
    return np.hypot(x - closest_x, y - closest_y)

def pctremain(x, exp_auc, both_auc):
    exp_ecdf = ecdf(exp_auc)
    both_ecdf = ecdf(both_auc)

    exp_remain = len(exp_auc) - exp_ecdf(x) * len(exp_auc)
    both_remain = len(both_auc) - both_ecdf(x) * len(both_auc)

    with np.errstate(divide='ignore', invalid='ignore'):
        final = np.divide(exp_remain, both_remain, out=np.zeros_like(exp_remain), where=both_remain != 0)

    return final


def cutoff_selection(exp_aucs, ctrl_aucs):
    both_aucs = np.concatenate((exp_aucs, ctrl_aucs))
    both_sorted_unique = np.sort(np.unique(both_aucs))

    # Vectorized pctremain computation
    pr_values = pctremain(both_sorted_unique, exp_aucs, both_aucs)
    pr_filtered = pr_values[pr_values < 1]

    if len(pr_filtered) == 0:
        raise ValueError("No values found where pctremain(x) < 1")

    max_pr = np.max(pr_filtered)
    x0_candidates = both_sorted_unique[pr_values == max_pr]

    if len(x0_candidates) == 0:
        raise ValueError("No x0 found where pctremain == max")

    x0 = x0_candidates[0]
    z = both_sorted_unique[both_sorted_unique <= x0]

    pr_z = pctremain(z, exp_aucs, both_aucs)
    midpoint = (pctremain(x0, exp_aucs, both_aucs) + np.min(pr_z)) / 2

    # Find z2: the z with minimum absolute diff from the midpoint
    z2 = z[np.argmin(np.abs(pr_z - midpoint))]

    if x0 != z2:
        z_post = z[z > z2]
        if len(z_post) == 0:
            z0 = z2
        else:
            half_dist = 0.5 * (np.max(z_post) - np.min(z_post))
            center = np.max(z_post) - half_dist
            z0 = z_post[np.argmin(np.abs(z_post - center))]
    else:
        z0 = x0

    # Compute empirical FDRs
    def fdr_threshold(x):
        return 1 - pctremain(x, exp_aucs, both_aucs)

    fdr_x0 = fdr_threshold(x0)
    fdr_z0 = fdr_threshold(z0)


    return x0, z0, fdr_x0, fdr_z0

def fast_ecdf(x, sorted_vals):
    return np.searchsorted(sorted_vals, x, side="right") / len(sorted_vals)

def load_raw_arrays(exp_array, ctrl_array):
    #Load in data into a numpy array.
    raw_exp = np.loadtxt(exp_array, delimiter="\t", dtype=[('auc', np.float32), ('length', np.uint32)])
    exp_data = np.column_stack((raw_exp['auc'], raw_exp['length'].astype(np.uint32)))
    
    raw_ctrl = np.loadtxt(ctrl_array, delimiter="\t", dtype=[('auc', np.float32), ('length', np.uint32)])
    ctrl_data = np.column_stack((raw_ctrl['auc'], raw_ctrl['length'].astype(np.uint32)))
    return exp_data, ctrl_data

########################################################################################################################


########################
# MOVING THINGS AROUND #
########################

def empirical_thresholding(exp, ctrl, norm):
    
    # Load in data as arrays
    exp_data, ctrl_data = load_raw_arrays(exp, ctrl)

    if norm:
        start_time = time.time()
        def base_function(data):
            # Sort input once
            sorted_input = data[data[:, 0].argsort()]
            sorted_aucs = sorted_input[:, 0]
            sorted_lengths = sorted_input[:, 1]
            
            # Calculate scaled AUCs and diff in one step
            scaled_aucs = sorted_aucs / np.max(sorted_aucs)
            line = np.linspace(1, 0, len(sorted_aucs))
            diff = line - scaled_aucs
            max_diff = np.max(diff)

            full_array = np.column_stack((sorted_aucs, sorted_lengths, scaled_aucs, line, diff))
            top_auc_array = full_array[full_array[:, 4] > 0.9 * max_diff]

            # Use vectorized operations to calculate distances
            distances = np.hypot(top_auc_array[:, 3] - top_auc_array[:, 2], top_auc_array[:, 0] - top_auc_array[:, 1])

            final_array = np.column_stack((top_auc_array, distances))
            return full_array, final_array

        def dist_threshold(array, cut_array):
            max_dist_idx = np.argmax(cut_array[:, 5])
            auc_at_max_dist = cut_array[max_dist_idx, 1]
            ninetieth_pct_idx = np.percentile(array[:, 0], 90) 

            return auc_at_max_dist if auc_at_max_dist > ninetieth_pct_idx else ninetieth_pct_idx
        
        # Base function and distance threshold for exp and ctrl
        exp_non_cut_arr, exp_thresh_arr = base_function(exp_data)
        ctrl_non_cut_arr, ctrl_thresh_arr = base_function(ctrl_data)
        
        exp_threshold_raw = dist_threshold(exp_non_cut_arr, exp_thresh_arr)
        ctrl_threshold_raw = dist_threshold(ctrl_non_cut_arr, ctrl_thresh_arr)

        def scale_ctrl_to_exp(array, threshold):
            # Filter values falling below the threshold
            subset = array[array <= threshold]
            if len(subset) < 2:
                raise ValueError("Not enough values under the threshold to estimate density")

            kde = gaussian_kde(subset, bw_method='silverman')
            xs = np.linspace(np.min(subset), np.max(subset), 500)
            ys = kde(xs)

            return xs[np.argmax(ys)]
        
        # Scale and calculate thresholds
        exp_mode = scale_ctrl_to_exp(exp_non_cut_arr[:, 0], exp_threshold_raw)
        ctrl_mode = scale_ctrl_to_exp(ctrl_non_cut_arr[:, 0], ctrl_threshold_raw)
        scaling_const = exp_mode / ctrl_mode
        ctrl_non_cut_arr[:, 0] *= scaling_const
        
        final_cutoff_exp, final_cutoff_ctrl, final_fdr_exp, final_fdr_ctrl = cutoff_selection(exp_non_cut_arr[:, 0], ctrl_non_cut_arr[:, 0])
        threshold = (final_cutoff_exp, final_cutoff_ctrl)
        fdr = (final_fdr_exp, final_fdr_ctrl)

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time} seconds")


        return threshold, fdr, scaling_const

    else:
        start_time = time.time()
        final_cutoff_exp, final_cutoff_ctrl, final_fdr_exp, final_fdr_ctrl = cutoff_selection(exp_data[:, 0], ctrl_data[:, 0])
        threshold = (final_cutoff_exp, final_cutoff_ctrl)
        fdr = (final_fdr_exp, final_fdr_ctrl)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time} seconds")
        return threshold, fdr

#TODO: Complete
def percentile_thresholding(exp, ctrl):
    threshold = "PLACEHOLDER"
    fdr = "PLACEHOLDER"


def threshold_calc(exp, ctrl, mode, normalize):
    print(f"threshold_calc started. Arguments are: {exp}, {ctrl}, {normalize}")
    if(mode == True):
        if(normalize == 'yes'):
            print("Empirical Thresholding with Normalization being performed.")
            threshold,fdr,norm_val = empirical_thresholding(exp, ctrl, True)
            norm = 1
            return threshold, fdr, norm, norm_val
        elif(normalize=='no'):
            print("Empirical Thresholding WITHOUT Normalization being performed.")
            threshold, fdr = empirical_thresholding(exp, ctrl, False)
            norm = 0
            return threshold, fdr, norm
        else:
            print("Selecting thresholding type and normalization has gone wrong")
    elif(mode == False):
        print("Percentile Thresholding being performed.")
        threshold,fdr = percentile_thresholding(exp, ctrl)
        norm = 0
        return threshold, fdr, norm

def check_ctrl(input_file):
    """Check if the input is a number or a file. If it's a file, load it."""
    global control_input
    
    try:
        # Try converting to float, if it works it's a number
        float(input_file)
        control_input = args.ctrl
        return False  # It's a number
    except ValueError:
        # If it fails, it's a file path
        file_path = Path(input_file).resolve()
        print(file_path)
        args.ctrl = file_path
        
        if file_path.exists():  # Check if the file exists
            control_input = np.loadtxt(file_path, delimiter="\t", dtype=[('col0', np.float32), ('col1', np.uint32)])
            print("File loaded.")
            return True  # It's a file and has been loaded
        else:
            exit()

#TODO: Complete
def write_output(output_prefix, cutoffs, fdrs, norm_type, norm_val=None):
    if args.norm == "yes":
        #Round thresholds
        rcutoffs = tuple(round(i, 4) for i in cutoffs)
        #Threshold file
        with open(f"{output_prefix}.threshold.txt", "w") as f:
            f.write(f"{rcutoffs[0]}\n{rcutoffs[1]}\n{norm_type}\n")
        #FDR file
        with open(f"{output_prefix}.fdr.txt", "w") as f:
            f.write(f"{fdrs[0]}\n{fdrs[1]}\n")
        #Normalization file
        with open(f"{output_prefix}.norm.txt", "w") as f:
            f.write(f"{norm_val}")
    else:
        #Round thresholds
        rcutoffs = tuple(round(i, 4) for i in cutoffs)
        #Threshold file
        with open(f"{output_prefix}.threshold.txt", "w") as f:
            f.write(f"{rcutoffs[0]}\n{rcutoffs[1]}\n{norm_type}\n")
        #FDR file
        with open(f"{output_prefix}.fdr.txt", "w") as f:
            f.write(f"{fdrs[0]}\n{fdrs[1]}\n")

########################################################################################################################

def main():
    ctrl_status = check_ctrl(args.ctrl)    
    if args.norm == "yes":
        threshold,fdr,norm_type, norm_val = threshold_calc(args.exp, args.ctrl, ctrl_status, args.norm)
        write_output(args.output,threshold, fdr, norm_type, norm_val)
    if args.norm == "no":
        threshold,fdr, norm_type = threshold_calc(args.exp, args.ctrl, ctrl_status, args.norm)
        write_output(args.output, threshold, fdr, norm_type)


main()
