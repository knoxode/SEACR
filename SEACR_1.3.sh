#!/usr/bin/env bash

#Treat unset variables as an error and exit immediately. Exit immediately  if any command exits with a non-zero exit status.
set -ue

#Checks if fewer than 5 arguments are provided. If true, gives help text.
# Checks if fewer than 5 arguments are provided. If true, gives help text.
if [ $# -lt 5 ]; then
	cat <<EOF
SEACR: Sparse Enrichment Analysis for CUT&RUN

Usage: bash SEACR_1.3.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDR threshold>] ["norm" | "non"] ["relaxed" | "stringent"] output prefix

Description of input fields:

Field 1: Target data bedgraph file in UCSC bedgraph format (https://genome.ucsc.edu/goldenpath/help/bedgraph.html) that omits regions containing 0 signal.

Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling. Alternatively, a numeric threshold n between 0 and 1 returns the top n fraction of peaks based on total signal within peaks.

Field 3: "norm" denotes normalization of control to target data, "non" skips this behavior. "norm" is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).

Field 4: "relaxed" uses a total signal threshold between the knee and peak of the total signal curve, and corresponds to the "relaxed" mode described in the text, whereas "stringent" uses the peak of the curve, and corresponds to "stringent" mode.

Field 5: Output prefix

Output file:
<output prefix>.auc.threshold.merge.bed (Bed file of enriched regions)

Output data structure:

<chr>	<start>	<end>	<AUC>	<max signal>	<max signal region>

Description of output fields:
Field 1: Chromosome
Field 2: Start coordinate
Field 3: End coordinate
Field 4: Total signal contained within denoted coordinates
Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates
Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal

Examples:
bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output
Calls enriched regions in target data using normalized IgG control track with stringent threshold

bash SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output
Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold

bash SEACR_1.3.sh target.bedgraph 0.01 non stringent output
Calls enriched regions in target data by selecting the top 1% of regions by area under the curve (AUC)
EOF
	exit 1
fi


#Generates random strings which will be used as placeholder names for temporary files.

random_string=$(head /dev/urandom | LC_CTYPE=C tr -dc A-Za-z0-9 | head -c 13; echo '')
random_string2=$(head /dev/urandom | LC_CTYPE=C tr -dc A-Za-z0-9 | head -c 13; echo '')

#Allows for use of easier-to-remember names in place of $1,$2, etc from bash arguments
exp=$1
ctrl=$2
norm=$3
height=$4

#Check if user has opted to forgo a control file.
if [[ $ctrl =~ ^[0-9]?+([.][0-9]+)?$ ]] || [[ $2 =~ ^[0-9]([.][0-9]+) ]] || [[ $2 =~ ^([.][0-9]+) ]]
then
	echo "Calling enriched regions without control file"
elif [[ -f $ctrl ]]
then
	echo "Calling enriched regions with control file"
else
	echo "$ctrl is not a number or a file"
	exit 1
fi

#Checks the normalization option provided by user.
if [[ $norm == "norm" ]]
then
	echo "Normalizing control to experimental bedgraph"
elif [[ $norm == "non" ]]
	then
	echo "Proceeding without normalization of control to experimental bedgraph"
else
	echo "Must specify \"norm\" for normalized or \"non\" for non-normalized data processing in third input"
	exit 1
fi

#Checks stringency option chosen by user.
if [[ $height == "relaxed" ]]
then
	echo "Using relaxed threshold"
elif [[ $height == "stringent" ]]
	then
	echo "Using stringent threshold"
else
	echo "Must specify \"stringent\" or \"relaxed\" in fourth input"
	exit 1
fi

echo "Creating experimental AUC file: $(date)"


generate_auc(){
  awk '
  # Initialize the state variable
  BEGIN { 
      s = 1 
    }

  # Process each line of input
  {
      #First line is header. Ignore, and move to line 2. Increment counter. 
      if (s == 1) {
          s++
      }
      #Assumed first real line of data.
      else if (s == 2) {
          #Checks signal is greater than zero.
          if ($4 > 0) {
              chr = $1         # Chromosome
              start = $2       # Start position
              stop = $3        # End position
              max = $4         # Max signal in this region
              coord = $1":"$2"-"$3  # Coordinate string in "chr:start-end" format
              auc = $4 * ($3 - $2)  # AUC for the region (signal * length of region)
              num = 1           # Counter for the number of regions in the block
              s++  # Move to the next state
          }
      }
      
      #Process lines 3 onward.
      else {
          if ($4 > 0) {
              #Check we are on the same chromosome, and are processing same region as previous line.
              if (chr == $1 && $2 == stop) {
                  num++  # Increment the region count
                  stop = $3  # Extend the end position of the block
                  auc += $4 * ($3 - $2)  # Add the signal * length for the new region to the AUC
                
                  # Update the maximum signal in this region if the current signal is higher
                  if ($4 > max) {
                      max = $4
                      coord = $1":"$2"-"$3  # Update the coordinate to reflect the region with the new max signal
                  } 

                  else if ($4 == max) {
                      split(coord, t, "-")  # Split the previous coord by the "-" to get start and end
                      coord = t[1] "-" $3  # Update the coordinate to reflect the new end position
                  }
              } 

              else {
                  # Print the results for the previous region block
                  print chr "\t" start "\t" stop "\t" auc "\t" max "\t" coord "\t" num
                  # Start a new block with the current region
                  chr = $1
                  start = $2
                  stop = $3
                  max = $4
                  coord = $1 ":" $2 "-" $3
                  auc = $4 * ($3 - $2)
                  num = 1
              }
          }
      }
  }'
  #Stream the output of awk into an auc.bed file
  "$1" > "$2".auc.bed
  #Stream the output of cut into a .auc file
  cut -f 4,7 "$2".auc.bed > "$2".auc
}

#Generate AUC for the sample
echo "Creating sample AUC file: $(date)"
generate_auc "$exp" "$random_string"

#If control_file is provided, generate control AUC.
if [[ -f $2 ]]
then
  echo "Creating control AUC file: $(date)"
  generate_auc "$ctrl" "$random_string2"
fi

#TODO: Remove after testing
path="dirname $0"
echo "$path"

auc_compare(){
	Rscript "$path"/SEACR_1.3.R --exp="$1".auc --ctrl="$2".auc --norm=yes --output="$1"
}

echo "Calculating optimal AUC threshold: $(date)"
path="dirname $0"
echo "$path"

if [[ -f $2 ]] && [[ $norm == "norm" ]]
then
	echo "Calculating threshold using normalized control: $(date)"
  auc_compare "$random_string" "$random_string2"
elif [[ -f $2 ]]
then
	echo "Calculating threshold using non-normalized control: $(date)"
  auc_compare "$random_string" "$random_string2"
else
	echo "Using user-provided threshold: $(date)"
  auc_compare "$random_string" "$ctrl"
fi


fdr="cat $random_string.fdr.txt | sed -n '1p'"
fdr2="cat $random_string.fdr.txt | sed -n '2p'"

#thresh=`cat $exp.threshold.txt`
thresh="cat $random_string.threshold.txt | sed -n '1p'"
thresh2="cat $random_string.threshold.txt | sed -n '2p'"
thresh3="cat $random_string.threshold.txt | sed -n '3p'"

echo "Creating thresholded feature file: $(date)"

if [[ $height == "relaxed" ]]
then
  echo "Empirical false discovery rate = $fdr2"
  awk -v value="$thresh2" -v value2="$thresh3" '$4 > value && $7 > value2 {print $0}' "$random_string".auc.bed | cut -f 1,2,3,4,5,6 > "$random_string".auc.threshold.bed
else
  echo "Empirical false discovery rate = $fdr"
  awk -v value="$thresh" -v value2="$thresh3" '$4 > value && $7 > value2 {print $0}' "$random_string".auc.bed | cut -f 1,2,3,4,5,6 > "$random_string".auc.threshold.bed
fi

if [[ -f $exp ]]
then
	if [[ $norm == "norm" ]] #If normalizing, multiply control bedgraph by normalization constant
	then
		constant="cat $random_string.norm.txt | sed -n '1p'"
		awk -v mult="$constant" 'BEGIN{OFS="\t"}; {$4=$4*mult; print $0}' "$random_string2".auc.bed | cut -f 1,2,3,4,5,6 > "$random_string2".auc2.bed
		mv "$random_string2".auc2.bed "$random_string2".auc.bed
	fi
	awk -v value="$thresh" '$4 > value {print $0}' "$random_string2".auc.bed > "$random_string2".auc.threshold.bed
fi

echo "Merging nearby features and eliminating control-enriched features: $(date)"


mean="awk '{s+=$3-$2; t++}END{print s/(t*10)}' $random_string.auc.threshold.bed"

if [[ -f $2 ]]; then
  awk -v value="$mean" '
    BEGIN { s = 1 }
    {
      if (s == 1) {
        chr = $1; start = $2; stop = $3;
        auc = $4; max = $5; coord = $6;
        s++;
      } else {
        if (chr == $1 && $2 < stop + value) {
          stop = $3;
          auc += $4;
          if ($5 > max) {
            max = $5;
            coord = $6;
          } else if ($5 == max) {
            split(coord, t, "-");
            split($6, u, "-");
            coord = t[1] "-" u[2];
          }
        } else {
          print chr "\t" start "\t" stop "\t" auc "\t" max "\t" coord;
          chr = $1; start = $2; stop = $3;
          auc = $4; max = $5; coord = $6;
        }
      }
    }
  ' "$random_string.auc.threshold.bed" | \
  bedtools intersect -wa -v \
    -a - \
    -b "$random_string2.auc.threshold.bed" \
    > "$5.auc.threshold.merge.bed"
else
  awk -v value="$mean" '
    BEGIN { s = 1 }
    {
      if (s == 1) {
        chr = $1; start = $2; stop = $3;
        auc = $4; max = $5; coord = $6;
        s++;
      } else {
        if (chr == $1 && $2 < stop + value) {
          stop = $3;
          auc += $4;
          if ($5 > max) {
            max = $5;
            coord = $6;
          } else if ($5 == max) {
            split(coord, t, "-");
            split($6, u, "-");
            coord = t[1] "-" u[2];
          }
        } else {
          print chr "\t" start "\t" stop "\t" auc "\t" max "\t" coord;
          chr = $1; start = $2; stop = $3;
          auc = $4; max = $5; coord = $6;
        }
      }
    }
  ' "$random_string.auc.threshold.bed" \
  > "$5.auc.threshold.merge.bed"
fi

if [[ $height == "relaxed" ]]
then
  cat "$5".auc.threshold.merge.bed > "$5".relaxed.bed
else
  cat "$5".auc.threshold.merge.bed > "$5".stringent.bed
fi

echo "Removing temporary files: $(date)"

rm "$random_string".auc.bed
rm "$random_string".auc
rm "$random_string".threshold.txt
rm "$random_string".auc.threshold.bed
rm "$random_string".fdr.txt
rm "$5".auc.threshold.merge.bed
if [[ -f $2 ]]
then
	rm "$random_string2".auc.bed
	rm "$random_string2".auc
	rm "$random_string2".auc.threshold.bed
fi
if [[ "$norm" == "norm" ]]
then
	rm -f "$random_string".norm.txt
fi
echo "Done: $(date)"
