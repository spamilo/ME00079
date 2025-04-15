#import necessary python packages
import pandas as pd
import numpy as np
import sys
from pathlib import Path
import os

#provide the path to local folder containing python functions for pbmc analysis
sys.path.insert(1, '/home/siina/python_scripts/moncyte_analysis/analysis_tools/pbmc_analysis_tools/')
#import functions for pbmc analysis
from DataPrep import DataPrep
from DataCleanup import DataCleanup

def cellp_to_single_cell(directory, filename, n_wells = 384, cutoff = 100, stain_a = None, stain_b=None,stain_c=None,stain_d=None):
    
    layouts_folder = os.path.join(directory,'layouts') # where your csv layout file is located
    results_folder = os.path.join(directory,'cellp_results') #where is your cellprofiler output folder located
    output_folder = os.path.join(directory,'out_results') #output data will be saved here 

    #extract experiment and run name from the filename
    parts = filename.split('_')
    exp_name = f'{parts[0]}_{parts[1]}'
    run_name = parts[0]

    fname = os.path.join(output_folder,f'single_cell_cleaned/{exp_name}_cell_clean.csv') # output filename
    log_file = os.path.join(directory, f'{fname[:-4]}_log.txt')

    with open(log_file, "w") as report_file:

        for f in [layouts_folder,results_folder]:
            dir_path = Path(f)
            if not dir_path.is_dir():
                raise FileNotFoundError(f"The directory '{f}' does not exist.")

        dir_path = Path(output_folder)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
            print(f"Directory '{output_folder}' created.")
            print(f"Directory '{output_folder}' created.", file=report_file)

        file_loc = os.path.join(results_folder,filename)
        print('Data location: ',file_loc)
        print('Data location: ',file_loc,file = report_file)

        print('Experiment name: ',exp_name)
        print('Experiment name: ',exp_name, file = report_file)
        
        print('Run name: ', run_name)
        print('Run name: ', run_name, file = report_file)

        experimentlayout = os.path.join(layouts_folder,f'{exp_name}_layout.csv')
        print('Plate layout: ',experimentlayout)
        print('Plate layout: ',experimentlayout, file = report_file)

        print(f'Number of wells = {n_wells}')
        print(f'Number of wells = {n_wells}', file = report_file)

        #initialise the class for transforming the dataframe to readable format
        my_frame = DataPrep(file_loc, exp_name)
        
        for i,stain in enumerate([stain_a,stain_b,stain_c,stain_d]):
            if stain:
                if i == 0:
                    my_frame.set_stain_a(stain)
                elif i == 1:
                    my_frame.set_stain_b(stain)
                elif i == 2:
                    my_frame.set_stain_c(stain)
                else:
                    my_frame.set_stain_d(stain)
        for s,n in zip([my_frame.stain_a, my_frame.stain_b, my_frame.stain_c, my_frame.stain_d], ['A','B','C','D']):
            print(f'Stain {n} = {s}')
            print(f'Stain {n} = {s}', file = report_file)

        #set the cut_off
        my_frame.set_cutoff(cutoff)
        print(f'Cut off = {cutoff}')
        print(f'Cut off = {cutoff}', file = report_file)
        dat = my_frame.org_calc(experimentlayout, n_wells)
        #remove outliers on a single-cell level
        clean_data = DataCleanup(5)
        cdat = clean_data.cl_fraction(dat,'cell_area')
        print('Removed ', dat.shape[0] - cdat.shape[0], ' cells out of ', dat.shape[0], ' based on cell area')
        print('Removed ', dat.shape[0] - cdat.shape[0], ' cells out of ', dat.shape[0], ' based on cell area', file = report_file)

        #cdat1 = clean_data.cl_well_fract(cdat, 'Size_DiIOrg')
        #print('removed ', cdat.shape[0] - cdat1.shape[0], ' cells out of ', cdat.shape[0], ' based on DiIOrg size')

        #save single-cell data with outliers removed 
        cdat.to_csv(fname)
        print('Cleaned single cell data saved to ' + fname)
        print('Cleaned single cell data saved to ' + fname, file = report_file)
    
    print(f"Logfile has been saved to {log_file}")

    return