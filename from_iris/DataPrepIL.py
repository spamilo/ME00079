# This set of functions is used for transforming the cell profiler output data into format readable
# by the rest of the scripts. The input file is the cellprofiler output csv file that ends with _Cells
# and layout file for the exxpetiment

import pandas as pd
import numpy as np
from pandas import DataFrame, Series

class DataPrep(object):
    # microscope resplution (pixel size) in micrometers, default is from OperaPhenix x40 objective
    pix_res = 0.149494
    # cutoff value based on which cells will be classified to monocytes of lymphocytes in square micrometers
    cutoff = 90
    # stains in different channels, a,b,c,d should be matched to your cellprofiler output
    stain_a = 'DAPI'
    stain_b = 'CellMask'
    stain_c = 'CD36'
    stain_d = 'DiI'
    # set of readouts that are not specific to a channel but will be included in the final output dataftame
    selection = ['Metadata_Filename','shortnames','experiment','well','Sample', 'treatment','time','Monocyte', 'um_cell_area']
    # info to be read from each cell of your layout file separated by space
    layout_opt = ['Sample', 'treatment']


    def __init__(self, path, exp_name):
        self.frame = pd.read_csv(path)
        self.exp_name = exp_name

    # the following functions are used to set the default values defined above to different ones

    def set_pixel_size(self, new_pix_res):
        self.pix_res = new_pix_res

    def set_cutoff(self, new_cutoff):
        self.cutoff = new_cutoff

    def set_selection(self, new_selection):
        self.selection = new_selection

    def set_stain_c(self, new_stain_c):
        self.stain_c = new_stain_c

    def set_stain_d(self, new_stain_d):
        self.stain_d = new_stain_d

    def set_layout_opt(self, new_layout_opt):
        self.layout_opt = new_layout_opt


    def three_eight_four_layout(self, layout_file):
        '''Reading patient and treatment info from a .csv file for an experiment in a 384-well plate'''
        template_384 = pd.read_csv('/home/ilmoncyte/Documents/analyses/cell_data_analyses/pbmc_analysis_tools/384wellplate.csv') # path to be updated on a different computer
        #print(template_384)
        experimentlayout = pd.read_csv(layout_file)
        #print(experimentlayout)
        platelay=pd.DataFrame({'well':template_384.values.ravel(), 'Names':experimentlayout.values.ravel()})
        print(platelay)
        platelay[self.layout_opt] = platelay['Names'].str.split(' ', expand=True)
        platelay=platelay.dropna()
        platelay=platelay.reset_index(drop=True)
        return(platelay)

    def nine_six_layout(self, layout_file):
        '''Reading patient and treatment info from a .csv file for an experiment in the 96-well plate'''
        template_96 = pd.read_csv('/home/pfistere/Documents/Laboratory/iryna_analysis_tools/pbmc_analysis_tools/96wellplate.csv')
        experimentlayout = pd.read_csv(layout_file)
        
        platelay=pd.DataFrame({'well':template_96.values.ravel(), 'Names':experimentlayout.values.ravel()})
        platelay[self.layout_opt] = platelay['Names'].str.split(' ', expand=True)
        platelay=platelay.dropna()
        platelay=platelay.reset_index(drop=True)
        return(platelay)

    def merge_from_csv(self, experimentlayout, pl):
        # pl variable is plate format, use 96 for 96-well format, anything else will use 384-well format
        ''' uses the platelayout dataframe generated in ..._layout functions, resolution, cutoff and stain names to make a readable dataframe with the cellinfo'''
        stain_a = self.stain_a
        stain_b = self.stain_b
        stain_c = self.stain_c
        stain_d = self.stain_d
        if pl==96:
            labels = self.nine_six_layout(experimentlayout)
        else:
            labels = self.three_eight_four_layout(experimentlayout)
        df = self.ren()
        exp_name = self.exp_name
        df['experiment'] = exp_name
        df['um_cell_area'] = df['cell_area']*(self.pix_res**2)
        df['Monocyte'] = np.where(df['um_cell_area']>self.cutoff, 'positive', 'negative')
        df['shortnames'] = df['Metadata_Filename'].str.slice(stop = 6)
        new_df = df.merge(labels, left_on = 'shortnames', right_on = 'well',how = 'outer')
        final_table_cols = self.selection + ['cell_area', stain_c+'Org', 'Size_'+stain_c+'Org', 'MeanInt'+stain_c+'_'+stain_c+'Org',
                            'MeanInt'+stain_d+'_'+stain_c+'Org','MeanMean'+stain_c+'_'+stain_c+'Org',
                            'MeanMean'+stain_d+'_'+stain_c+'Org', 'MeanMed'+stain_c+'_'+stain_c+'Org',
                            'MeanMed'+stain_d+'_'+stain_c+'Org', 'Mean'+stain_a, 'Mean'+stain_b, 'Mean'+stain_c,
                            'Mean'+stain_d, 'Int'+stain_a, 'Int'+stain_b, 'Int'+stain_c, 'Int'+stain_d, 'Med'+stain_a,
                            'Med'+stain_b,'Med'+stain_c, 'Med'+stain_d, 'Mean'+stain_c+'Nuclei', 'Med'+stain_c+'Nuclei',
                            stain_d+'Org', 'Size_'+stain_d+'Org', 'MeanInt'+stain_c+'_'+stain_d+'Org',
                            'MeanInt'+stain_d+'_'+stain_d+'Org', 'MeanMean'+stain_c+'_'+stain_d+'Org',
                            'MeanMean'+stain_d+'_'+stain_d+'Org', 'MeanMed'+stain_c+'_'+stain_d+'Org',
                            'MeanMed'+stain_d+'_'+stain_d+'Org']
        new_df = new_df[new_df.columns.intersection(final_table_cols)]
        return new_df

    def combine_oa_files(self, path, spot_file_name, cell_file_name):
        '''Makes a single csv file from two Object Analyzer output files. It has the same column names as cellprofiler cell measurement file'''
        self.path = path
        self.spot_file_name = spot_file_name
        self.cell_file_name = cell_file_name
        spot_df = pd.read_csv(self.path + self.spot_file_name)
        cell_df = pd.read_csv(self.path + self.cell_file_name)
        mean_spot = spot_df.groupby(['ImageName', 'Parent_MetaData_Parent::Cells']).mean()
        mean_spot = mean_spot.rename(columns = {'Area_AreaShape_1pxM1m':'Mean_SpotsC_AreaShape_Area',
                                                'MeanIntensity_IntensityGray_ch2::LD' : 'Mean_SpotsC_Intensity_MeanIntensity_C'})
        #print(mean_spot.head(3))
        count_spot = spot_df.groupby(['ImageName', 'Parent_MetaData_Parent::Cells']).Area_AreaShape_1pxM1m.count()
        #print(count_spot.head(3))
        comb = pd.concat([mean_spot, count_spot], axis = 1)
        comb = comb.rename(columns = {'Area_AreaShape_1pxM1m':'Children_SpotsC_Count'})
        #print(comb.head(3))
        cell_df = cell_df.rename(columns = {'Object index_MetaData_Parent::Cells':'Parent_MetaData_Parent::Cells'})
        cell_df = cell_df.set_index(['ImageName', 'Parent_MetaData_Parent::Cells'])
        comb2 = pd.concat([cell_df, comb], axis = 1).fillna(0).reset_index()
        comb2 = comb2.rename(columns = {'Area_AreaShape_1pxM1m': 'cell_area', 'Parent_MetaData_Parent::Cells':'Cell_ID',
                                            'ImageName':'Metadata_Filename'})
        comb2 = comb2.loc[comb2['cell_area'] != 0]
        comb2.to_csv(self.path+self.exp_name+'_OAcombined.csv')
        return comb2

    def org_calc(self, experimentlayout, pl):
        '''Calculates um size of organelles and calculates cellular org area and intensity based on merge info'''

        ren_df = self.merge_from_csv(experimentlayout, pl)
        stain_c = self.stain_c
        stain_d = self.stain_d

        if stain_c+'Org' in ren_df.columns:
            ren_df[stain_c+'Org_pos'] = ren_df[stain_c+'Org' ] > 0
            ren_df[stain_c+'Org_pos'] = ren_df[stain_c+'Org_pos'].replace({True: 'positive', False: 'negative'})

        if stain_d+'Org' in ren_df.columns:
            ren_df[stain_d+'Org_pos'] = ren_df[stain_d+'Org' ] > 0
            ren_df[stain_d+'Org_pos'] = ren_df[stain_d+'Org_pos'].replace({True: 'positive', False: 'negative'})

        if 'Size_'+stain_c+'Org' in ren_df.columns:
            ren_df['um_'+stain_c+'Org_size'] = ren_df['Size_'+stain_c+'Org']*(self.pix_res**2)
            if stain_c+'Org' in ren_df.columns:
                ren_df['total_'+stain_c+'Org_size'] = ren_df[stain_c+'Org']*ren_df['Size_'+stain_c+'Org']
                ren_df['um_total_'+stain_c+'Org_size'] = ren_df[stain_c+'Org']*ren_df['um_'+stain_c+'Org_size']

        if stain_c+'Org' in ren_df.columns:
            if 'MeanMean'+stain_c+'_'+stain_c+'Org' in ren_df.columns:
                ren_df['total_'+stain_c+'_meanmean_of_'+stain_c+'Org'] = ren_df[stain_c+'Org']*ren_df['MeanMean'+stain_c+'_'+stain_c+'Org']
            if 'MeanMed'+stain_c+'_'+stain_c+'Org' in ren_df.columns:
                ren_df['total_'+stain_c+'_meanmed_of_'+stain_c+'Org'] = ren_df[stain_c+'Org']*ren_df['MeanMed'+stain_c+'_'+stain_c+'Org']
            if 'MeanInt'+stain_c+'_'+stain_c+'Org' in ren_df.columns:
                ren_df['total_'+stain_c+'_meanint_of_'+stain_c+'Org'] = ren_df[stain_c+'Org']*ren_df['MeanInt'+stain_c+'_'+stain_c+'Org']
            if 'MeanInt'+stain_d+'_'+stain_c+'Org' in ren_df.columns:
                ren_df['total_'+stain_d+'_meanint_of_'+stain_c+'Org'] = ren_df[stain_c+'Org']*ren_df['MeanInt'+stain_d+'_'+stain_c+'Org']
            if 'MeanMean'+stain_d+'_'+stain_c+'Org' in ren_df.columns:
                ren_df['total_'+stain_d+'_meanmean_of_'+stain_c+'Org'] = ren_df[stain_c+'Org']*ren_df['MeanMean'+stain_d+'_'+stain_c+'Org']
            if 'MeanMed'+stain_d+'_'+stain_c+'Org' in ren_df.columns:
                ren_df['total_'+stain_d+'_meanmed_of_'+stain_c+'Org'] = ren_df[stain_c+'Org']*ren_df['MeanMed'+stain_d+'_'+stain_c+'Org']

        if 'Size_'+stain_d+'Org' in ren_df.columns:
            ren_df['um_'+stain_d+'Org_size'] = ren_df['Size_'+stain_d+'Org']*(self.pix_res**2)
            if stain_d+'Org'in ren_df.columns:
                ren_df['total_'+stain_d+'Org_size'] = ren_df[stain_d+'Org']*ren_df['Size_'+stain_d+'Org']
                ren_df['um_total_'+stain_d+'Org_size'] = ren_df[stain_d+'Org']*ren_df['um_'+stain_d+'Org_size']

        if stain_d+'Org'in ren_df.columns:
            if 'MeanInt'+stain_c+'_'+stain_d+'Org' in ren_df.columns:
                ren_df['total_'+stain_c+'_meanint_of_'+stain_d+'Org'] = ren_df[stain_d+'Org']*ren_df['MeanInt'+stain_c+'_'+stain_d+'Org']
            if 'MeanMean'+stain_c+'_'+stain_d+'Org' in ren_df.columns:
                ren_df['total_'+stain_c+'_meanmean_of_'+stain_d+'Org'] = ren_df[stain_d+'Org']*ren_df['MeanMean'+stain_c+'_'+stain_d+'Org']
            if 'MeanMed'+stain_c+'_'+stain_d+'Org' in ren_df.columns:
                ren_df['total_'+stain_c+'_meanmed_of_'+stain_d+'Org'] = ren_df[stain_d+'Org']*ren_df['MeanMed'+stain_c+'_'+stain_d+'Org']
            if 'MeanInt'+stain_d+'_'+stain_d+'Org' in ren_df.columns:
                ren_df['total_'+stain_d+'_meanint_of_'+stain_d+'Org'] = ren_df[stain_d+'Org']*ren_df['MeanInt'+stain_d+'_'+stain_d+'Org']
            if 'MeanMean'+stain_d+'_'+stain_d+'Org' in ren_df.columns:
                ren_df['total_'+stain_d+'_meanmean_of_'+stain_d+'Org'] = ren_df[stain_d+'Org']*ren_df['MeanMean'+stain_d+'_'+stain_d+'Org']
            if 'MeanMed'+stain_d+'_'+stain_d+'Org' in ren_df.columns:
                ren_df['total_'+stain_d+'_meanmed_of_'+stain_d+'Org'] = ren_df[stain_d+'Org']*ren_df['MeanMed'+stain_d+'_'+stain_d+'Org']

        return ren_df

    def ren(self):
        '''This function renames the cell profiler column names of the dataframe'''
        stain_a = self.stain_a
        stain_b = self.stain_b
        stain_c = self.stain_c
        stain_d = self.stain_d
        x = self.frame
        x = x.rename(columns= {'AreaShape_Area' : 'cell_area',
                             'Children_SpotsC_Count' : stain_c+'Org',
                             'Mean_SpotsC_AreaShape_Area' : 'Size_'+stain_c+'Org',
                             'Mean_SpotsC_Intensity_IntegratedIntensity_C' : 'MeanInt'+stain_c+'_'+stain_c+'Org',
                             'Mean_SpotsC_Intensity_IntegratedIntensity_D' : 'MeanInt'+stain_d+'_'+stain_c+'Org',
                             'Mean_SpotsC_Intensity_MeanIntensity_C' : 'MeanMean'+stain_c+'_'+stain_c+'Org',
                             'Mean_SpotsC_Intensity_MeanIntensity_D' : 'MeanMean'+stain_d+'_'+stain_c+'Org',
                             'Mean_SpotsC_Intensity_MedianIntensity_C': 'MeanMed'+stain_c+'_'+stain_c+'Org',
                             'Mean_SpotsC_Intensity_MedianIntensity_D': 'MeanMed'+stain_d+'_'+stain_c+'Org',
                             'Intensity_MeanIntensity_A' : 'Mean'+stain_a,
                             'Intensity_MeanIntensity_B' : 'Mean'+stain_b,
                             'Intensity_MeanIntensity_C' : 'Mean'+stain_c,
                             'Intensity_MeanIntensity_D' : 'Mean'+stain_d,
                             'Intensity_IntegratedIntensity_A': 'Int'+stain_a,
                             'Intensity_IntegratedIntensity_B': 'Int'+stain_b,
                             'Intensity_IntegratedIntensity_C': 'Int'+stain_c,
                             'Intensity_IntegratedIntensity_D': 'Int'+stain_d,
                             'Intensity_MedianIntensity_A': 'Med'+stain_a,
                             'Intensity_MedianIntensity_B': 'Med'+stain_b,
                             'Intensity_MedianIntensity_C': 'Med'+stain_c,
                             'Intensity_MedianIntensity_D': 'Med'+stain_d,
                             'Mean_Nuclei_Intensity_MeanIntensity_C': 'Mean'+stain_c+'Nuclei',
                             'Mean_Nuclei_Intensity_MedianIntensity_C': 'Med'+stain_c+'Nuclei',
                             'Children_SpotsD_Count': stain_d+'Org',
                             'Mean_SpotsD_AreaShape_Area' : 'Size_'+stain_d+'Org',
                             'Mean_SpotsD_Intensity_IntegratedIntensity_C' : 'MeanInt'+stain_c+'_'+stain_d+'Org',
                             'Mean_SpotsD_Intensity_IntegratedIntensity_D' : 'MeanInt'+stain_d+'_'+stain_d+'Org',
                             'Mean_SpotsD_Intensity_MeanIntensity_C' : 'MeanMean'+stain_c+'_'+stain_d+'Org',
                             'Mean_SpotsD_Intensity_MeanIntensity_D' : 'MeanMean'+stain_d+'_'+stain_d+'Org',
                             'Mean_SpotsD_Intensity_MedianIntensity_C': 'MeanMed'+stain_c+'_'+stain_d+'Org',
                             'Mean_SpotsD_Intensity_MedianIntensity_D': 'MeanMed'+stain_d+'_'+stain_d+'Org',
                             })
        return x
