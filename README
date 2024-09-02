# Automated Smear Analysis for Fragment Analyzer Data

## Functionality

This script is designed to automatically define the edges of the target peak to calculate the percent RNA integrity from electropherogram data exported from ProSize.
The script is able to handle multiple sample runs at a time by placing this folder `AutoSmear/` in the same directory as the folders containing the output from each FA run (as diagrammed in `Usage`)

## Usage

1. Upload the electropherogram data (.raw) file obtained from the fragment analyzer to ProSize, then export the data to .csv format as a single file with "size" as the x-axis. 

2. Place this folder `AutoSmear/` in a directory containing one or more of the folders obtained from exporting from ProSize as described in step 1. The directory layout is described in the figure below.

            .
            └── folder/
                ├── **AutoSmear/**
                │   ├── AutoSmear.py
                │   └── AutoSmear.sh
                ├── ProSize_run1/
                │   └── Electropherogram.csv
                ├── ProSize_run2/
                │   └── Electropherogram.csv
                └── ProSize_run3/
                    └── Electropherogram.csv


3. Open a terminal on Mac or another command line interface. With the current directory `AutoSmear/` set as the working directory, run the `AutoSmear.sh` script with the following command:

        bash AutoSmear.sh

4. The results of the script are written to a .csv file as a table containing ID, left gate, right gate, and percent integrity. The csv file can be found in a folder called `smear_results/`, which is created by the script

##Outputs

The outputs of the script are written into a folder `smear_results/`, which can be found in the same directory as this folder `AutoSmear/`

There are two outputs from this script:

    A .csv file containing the smear analysis results formatted as seen below

        | Lane | ID        | Left_Gate | Right_Gate | Percent_Total |
        |------|-----------|-----------|------------|---------------|
        | A1   | sample1_1 | 9000      | 11000      | 80.00         |
        | A2   | sample1_2 | 9000      | 11000      | 70.00         |
        | A3   | sample1_3 | 9000      | 11000      | 60.00         |

    A .jpg file containing plots for each capillary, denoting where the peak boundaries were set for integrity calculations.
    This is essential for visually confirming proper peak identification and identifying poor quality capillary runs.



