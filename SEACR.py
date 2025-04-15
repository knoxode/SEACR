import argparse
from pathlib import Path
import numpy as np


#   data<-read.table(datasource)
#   datavec<-data$V1
#   datamax<-data$V2
#
#   dist2d<-function(a,b,c){v1<- b - c; v2<- a - b; m<-cbind(v1,v2); d<-det(m)/sqrt(sum(v1*v1))}
#
#   #Transform data into a 0 to 1 number space.
#   dataframe<-data.frame(count=seq(1,0,length=length(datavec)), quant=sort(datavec,decreasing=TRUE)/max(datavec), value=sort(datavec,decreasing=TRUE))
#
#   dataframe$diff<-abs(dataframe$count - dataframe$quant)
#   dataframe<-dataframe[dataframe$diff > 0.9*max(dataframe$diff),]
#   dataframe$dist<-apply(dataframe,1,function(x) dist2d(c(x[1],x[2]),0,1))
#
#   return(list(dataframe=dataframe, datavec=datavec, datamax=datamax))
#

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

#RANDOM CRAP THAT MAY BE OF USE LATER


# max_sig_exp = np.max(experiment_input[:,0])
# max_sig_ctrl = np.max(control_input[:,0])




########################


def dist2d(a, b, c):
    """Perpendicular distance from point a to line bc."""
    v1 = np.subtract(b, c)
    v2 = np.subtract(a, b)
    d = np.linalg.det([v1, v2]) / np.linalg.norm(v1)
    return abs(d)



#TODO: Complete
def empirical_thresholding(exp, ctrl, norm):
    #Load in data into a numpy array.
    exp_data = np.loadtxt(exp, delimiter="\t", dtype=[('auc', np.float32), ('length', np.uint32)])
    ctrl_data = np.loadtxt(norm, delimiter="\t", dtype=[('auc', np.float32), ('length', np.uint32)])

    if (norm == True):

        def base_function(data):

            #Find max auc from each region bin
            max_auc = np.max(data[:,0])
            #Scale AUCs between 0 and 1
            scaled_aucs = data[:, 0] / np.max(data[:, 0])
            #Sorts the lengths of bin regions
            sorted_lengths = np.sort(data[:,1], kind='quicksort')
            #Compares the length of signal blocks to the max AUC
            fraction_of_max = sorted_lengths/max_auc
            #Finds the difference between scaled aucs and fraction of max for each element
            difference = np.max(scaled_aucs - fraction_of_max)
            #Rejoining elements to create an array
            full_array = np.column_stack((scaled_aucs, fraction_of_max, sorted_lengths))
    
            top_auc_array = full_array[full_array[:, 2] >= 0.9 * difference]

            
            # Initialize an empty array to store distances
            distances = np.zeros(data.shape[0])
            # Calculate distance for each row
            for i in range(data.shape[0]):
                distances[i] = dist2d(data[i, 0], data[i, 1], 0)
            # Add the 'dist' column to final array
            final_array = np.column_stack((top_auc_array, distances))

            # Show the result
            return final_array

        exp_curve = base_function(exp_data)
        ctrl_curve = base_function(ctrl_data)

        #TODO: Calculate 

        def dist_threshold(array):
            max_dist_idx = np.argmax(array[:, 3])
            auc_at_max_dist = array[max_dist_idx,1]
            ninetieth_pct_idx = np.percentile(array[:,3], 90) 

            if(auc_at_max_dist > ninetieth_pct_idx):
                value = auc_at_max_dist
                return value
            else:
                value = ninetieth_pct_idx
                return value
        
        dist_threshold(exp_curve)
        dist_threshold(ctrl_curve)

        threshold = "PLACEHOLDER"
        fdr = "PLACEHOLDER"
        write_output(fdr, threshold) 


    elif (norm == False):
        threshold = "PLACEHOLDER"
        fdr = "PLACEHOLDER"
        write_output(fdr, threshold)


    else:
        print("Normalization parameter has not been set to either \"Yes\" or \"No\". Please set this correctly, and try again.")



#TODO: Complete
def percentile_thresholding(exp, ctrl):
    threshold = "PLACEHOLDER"
    fdr = "PLACEHOLDER"
    write_output(fdr, threshold)


def threshold_calc(exp, ctrl, mode, normalize):
    if(mode == "empirical"):
        if(normalize == "yes"):
            empirical_thresholding(exp, ctrl, True)
        elif(normalize=='no'):
            empirical_thresholding(exp, ctrl, False)
        else:
            print("Selecting thresholding type and normalization has gone wrong")
    elif(mode == "percentile"):
        percentile_thresholding(exp, ctrl)

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
        file_path = Path(input_file)
        
        if file_path.exists():  # Check if the file exists
            control_input = np.loadtxt(file_path, delimiter="\t", dtype=[('col0', np.float32), ('col1', np.uint32)])
            return True  # It's a file and has been loaded
        else:
            print(f"Error: The control file {input_file} does not exist.")
            exit()

#TODO: Complete
def write_output(threshold, fdr):
    #write thresholds to file using args.output prefix
    #write fdrs to file using args.output prefix
    pass

def main():
    ctrl_status = check_ctrl(args.ctrl)    
    threshold_calc(args.exp, args.ctrl, ctrl_status, args.norm)


main()
