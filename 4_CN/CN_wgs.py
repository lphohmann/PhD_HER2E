# !/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Title:     Copy number alterations in the HER2E subtype (WGS based)
Date:      01.08.2022
Author:    Lennart Hohmann

Usage:     CN_wgs.py ASCAT_file1 ASCAT_file2 ... output_file

TODO:

'''

##############################################################################
#%%

''' import of modules '''

#import argparse
import pandas as pd
from functools import reduce

# wdir
wdir = "~/Desktop/MTP_project/"

##############################################################################
#%%

''' user input : Ill do that stuff later, for now i will just hardcode the file loading'''

'''
parser = argparse.ArgumentParser(description="This program takes a fasta file containing one prokaryotic genome and predicts genes based on the presence of a Shine-Dargarno sequence before an open reading frame (ORF).")
parser.add_argument('-g', '--genomefile', help="The file containing the genome sequence in fasta format.", required=True)
parser.add_argument('-o', '--outputfile', help="The desired name of the output file containing the predicted genes. (default=geneprediction_output)", default='geneprediction_output')
parser.add_argument('-ml', '--minORFlength', help="The minimum ORF length of predicted genes. (default=300)", default=300, type = int)
parser.add_argument('-sc', '--startcodons', help="The possible start codons. Input like: -sc ATG GTG TTG ; (default=ATG)", nargs='+', default=["ATG"], type=str)
parser.add_argument('-ec', '--stopcodons', help="The possible stop codons. Input like: -ec TAA TAG TGA ; (default=TAA TAG TGA)", nargs='+', default=["TAA","TAG","TGA"], type=str)
parser.add_argument('-sd', '--shinedalgarnosequence', help="The Shine-Dalgarno sequence. (default=AGGAGG)", default="AGGAGG", type=str)
parser.add_argument('-me', '--maxalignerrors', help="The maximal number of errors that are allowed in the pairwise alignemnt when determining whether a Shine-Dalgarno sequence is present or not. (default=1)", default=1, type=int)
parser.add_argument('-f', '--fastaoutput', help="Output a file with the DNA sequences of the predicted genes in fasta format", action='store_true')

args = parser.parse_args()
args
# saving user input in specific variables
start_codons = set([codon.upper() for codon in args.startcodons]) # convert to uppercase and set
stop_codons = set([codon.upper() for codon in args.stopcodons]) # convert to uppercase and set
input_file = args.genomefile
output_file = args.outputfile
min_length = args.minORFlength
SD_seq = args.shinedalgarnosequence.upper()
allowed_errors = args.maxalignerrors
fasta_output = args.fastaoutput
'''

###############################################################################

#%%

''' load all files as panda dataframes: i hardcode everything for now; will change later with os etc.'''

# load readme for colnames
readme_file = pd.read_csv(wdir+'Data/SCAN_B/WGS/ASCATreadme.tsv', sep='\t', header=None)
#readme_file.info()
colnames = list(readme_file.iloc[1:9,0].str.strip())


# load ascat files

S000763_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S000763_l_d_a.copynumber.caveman.csv",header=None)
S000763_l_d_a.columns = colnames

S001143_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S001143_l_d_a.copynumber.caveman.csv",header=None)
S001143_l_d_a.columns = colnames 

S001347_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S001347_l_d_a.copynumber.caveman.csv",header=None)
S001347_l_d_a.columns = colnames

S002097_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S002097_l_d_a.copynumber.caveman.csv",header=None)
S002097_l_d_a.columns = colnames

S002236_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S002236_l_d_a.copynumber.caveman.csv",header=None)
S002236_l_d_a.columns = colnames

S002369_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S002369_l_d_a.copynumber.caveman.csv",header=None)
S002369_l_d_a.columns = colnames 

S002605_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S002605_l_d_a.copynumber.caveman.csv",header=None)
S002605_l_d_a.columns = colnames

S003100_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S003100_l_d_a.copynumber.caveman.csv",header=None)
S003100_l_d_a.columns = colnames

S003516_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S003516_l_d_a.copynumber.caveman.csv",header=None)
S003516_l_d_a.columns = colnames

S004149_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S004149_l_d_a.copynumber.caveman.csv",header=None)
S004149_l_d_a.columns = colnames

S004725_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S004725_l_d_a.copynumber.caveman.csv",header=None)
S004725_l_d_a.columns = colnames

S005529_l_d_a = pd.read_csv(wdir+"Data/SCAN_B/WGS/ASCAT/S005529_l_d_a.copynumber.caveman.csv",header=None)
S005529_l_d_a.columns = colnames


ascatfile_list = [S000763_l_d_a, 
                  S001143_l_d_a, 
                  S001347_l_d_a, 
                  S002097_l_d_a,
                  S002236_l_d_a, 
                  S002369_l_d_a, 
                  S002605_l_d_a, 
                  S003100_l_d_a, 
                  S003516_l_d_a, 
                  S004149_l_d_a, 
                  S004725_l_d_a, 
                  S005529_l_d_a]

ascatname_list = ["S000763_l_d_a", 
                  "S001143_l_d_a", 
                  "S001347_l_d_a", 
                  "S002097_l_d_a",
                  "S002236_l_d_a", 
                  "S002369_l_d_a", 
                  "S002605_l_d_a", 
                  "S003100_l_d_a", 
                  "S003516_l_d_a", 
                  "S004149_l_d_a", 
                  "S004725_l_d_a", 
                  "S005529_l_d_a"]


# load the probe file
probe_file = pd.read_csv(wdir+'Data/SCAN_B/WGS/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv', sep='\t', header=0)
probe_file = probe_file.iloc[:,0:3]
probe_file.columns = ["ProbeID","Chr","Position"]
#probe_file.info()
probe_file['Chr'] = probe_file['Chr'].map(lambda x: x.lstrip('chr'))
#probe_file["Chr"] = pd.to_numeric(probe_file["Chr"]) #doesnt work due to x CHR
#probe_file

# load file with sample ploidy info
ploidy_file = pd.read_excel(wdir+"Data/SCAN_B/WGS/HER2enrichedSummary.xlsx",header=0,sheet_name='SampleSummary')

ploidy_file = ploidy_file[["Tumour", "Ploidy"]]

###############################################################################

#%%

''' defining functions '''

# function that creates a CN matrix (probeID | Chr | Position | CNstate_sample) for a single chromosome with a method defining what type of matrix is created (total, major, minor)
def make_chrmatrix(chr_probes, chr_segments, matrix_type, tumploidy):
    # here it has ProbeID, Chr, Position
    chromatrix = chr_probes
    # create new column
    chromatrix['CN_state'] = '' 
    
    
    # types: major; minor; total; gainloss; amplification; LOH
    if matrix_type == "major":
        # look up the segments now
        for seg in chr_segments['Incremental count of segment']:
            # segment data
            seg_data = chr_segments.loc[chr_segments['Incremental count of segment']==seg]
            seg_end = seg_data.iloc[0]["stop position"]
            seg_start = seg_data.iloc[0]["start position"]
            seg_CNstate = seg_data.iloc[0]['Total CN tumour'] - seg_data.iloc[0]['Minor CN tumour']
            # map probes to segments now
            chromatrix.loc[(chromatrix['Position'] >= seg_start) & (chromatrix['Position'] <= seg_end), 'CN_state'] = seg_CNstate
            
            
    elif matrix_type == "minor":
        # look up the segments now
        for seg in chr_segments['Incremental count of segment']:
            # segment data
            seg_data = chr_segments.loc[chr_segments['Incremental count of segment']==seg]
            seg_end = seg_data.iloc[0]["stop position"]
            seg_start = seg_data.iloc[0]["start position"]
            seg_CNstate = seg_data.iloc[0]['Minor CN tumour']
            # map probes to segments now
            chromatrix.loc[(chromatrix['Position'] >= seg_start) & (chromatrix['Position'] <= seg_end), 'CN_state'] = seg_CNstate
            
            
    elif matrix_type == "total":
        # look up the segments now
        for seg in chr_segments['Incremental count of segment']:
            # segment data
            seg_data = chr_segments.loc[chr_segments['Incremental count of segment']==seg]
            seg_end = seg_data.iloc[0]["stop position"]
            seg_start = seg_data.iloc[0]["start position"]
            seg_CNstate = seg_data.iloc[0]['Total CN tumour']
            # map probes to segments now
            chromatrix.loc[(chromatrix['Position'] >= seg_start) & (chromatrix['Position'] <= seg_end), 'CN_state'] = seg_CNstate
            
            
    elif matrix_type == "gainloss":
        # look up the segments now
        for seg in chr_segments['Incremental count of segment']:
            # segment data
            seg_data = chr_segments.loc[chr_segments['Incremental count of segment']==seg]
            seg_end = seg_data.iloc[0]["stop position"]
            seg_start = seg_data.iloc[0]["start position"]
            if seg_data.iloc[0]['Total CN tumour'] > (tumploidy + 0.6): # gain
                seg_CNstate = 1
            elif seg_data.iloc[0]['Total CN tumour'] < (tumploidy - 0.6): # loss
                seg_CNstate = -1
            else: seg_CNstate = 0 # neutral
            # map probes to segments now
            chromatrix.loc[(chromatrix['Position'] >= seg_start) & (chromatrix['Position'] <= seg_end), 'CN_state'] = seg_CNstate    
            
            
    elif matrix_type == "amplification":
        # look up the segments now
        for seg in chr_segments['Incremental count of segment']:
            # segment data
            seg_data = chr_segments.loc[chr_segments['Incremental count of segment']==seg]
            seg_end = seg_data.iloc[0]["stop position"]
            seg_start = seg_data.iloc[0]["start position"]
            if seg_data.iloc[0]['Total CN tumour'] > (tumploidy*4): # focal amp
                seg_CNstate = 2
            elif seg_data.iloc[0]['Total CN tumour'] > (tumploidy*2): # amp
                seg_CNstate = 1
            else: seg_CNstate = 0 # no amp
            # map probes to segments now
            chromatrix.loc[(chromatrix['Position'] >= seg_start) & (chromatrix['Position'] <= seg_end), 'CN_state'] = seg_CNstate
            
            
    elif matrix_type == "LOH":
        # look up the segments now
        for seg in chr_segments['Incremental count of segment']:
            # segment data
            seg_data = chr_segments.loc[chr_segments['Incremental count of segment']==seg]
            seg_end = seg_data.iloc[0]["stop position"]
            seg_start = seg_data.iloc[0]["start position"]
            if seg_data.iloc[0]['Minor CN tumour']==0: # LOH
                seg_CNstate = 1
            else: seg_CNstate = 0 # not LOH
            # map probes to segments now
            chromatrix.loc[(chromatrix['Position'] >= seg_start) & (chromatrix['Position'] <= seg_end), 'CN_state'] = seg_CNstate
    else: 
        print("check matrix_type input")
    return chromatrix


#%% 
'''
# function that creates a matrix (probeID | Chr | Position | CNstate_sample) for a single chromosome
def make_chromatrix(chr_probes, chr_segments):
    # here it has ProbeID, Chr, Position
    chromatrix = chr_probes
    # create new column
    chromatrix['CN_state'] = '' 
    # look up the segments now
    for seg in chr_segments['Incremental count of segment']:
        # segment data
        seg_data = chr_segments.loc[chr_segments['Incremental count of segment']==seg]
        seg_end = seg_data.iloc[0]["stop position"]
        seg_start = seg_data.iloc[0]["start position"]
        seg_CNstate = seg_data.iloc[0]['Total CN tumour'] - seg_data.iloc[0]['Minor CN tumour']
        # map probes to segments now
        chromatrix.loc[(chromatrix['Position'] >= seg_start) & (chromatrix['Position'] <= seg_end), 'CN_state'] = seg_CNstate
    return chromatrix
'''

# function that that creates a matrix (probeID | Chr | Position | CNstate_sample) for a single sample
def make_sampmatrix(probe_file, ascat_file, sampleID, matrix_type, tumploidy):
    chr_set = set(ascat_file.iloc[:,1])
    res_list = []
    for chr in chr_set:
        # relevant probes
        chr_probes = probe_file.loc[probe_file['Chr'] == chr]
        # relevant segments
        chr_segments = ascat_file.loc[ascat_file['Chr'] == chr]
        # get chr matrix
        chromatrix = make_chrmatrix(chr_probes, chr_segments,matrix_type,tumploidy)

        # save in list
        res_list.append(chromatrix)
    
    # try to make one dataframe out of the list
    sampmatrix = pd.concat(res_list)
    # correct data type
    #sampmatrix[['Position', 'CN_state']] = sampmatrix[['Position', 'CN_state']].apply(pd.to_numeric)
    # HERE REANME CN_STATE ACCORDING TO ASCAT FILE NAME
    sampmatrix.rename(columns={'CN_state': sampleID}, inplace=True)

    return sampmatrix

# function that that creates a matrix (probeID | Chr | Position | CNstate_sample) for a all samples
def make_finmatrix(probe_file, ascatfile_list, ascatname_list, ploidy_file ,matrix_type):
    fin_list = []
    i = 0
    for ascat_file in ascatfile_list:
        # get tumor ploidy for that sample
        samp_index = ploidy_file.index[ploidy_file['Tumour'] == ascatname_list[i]].tolist()[0] # get row index for sample
        tumploidy = ploidy_file.iloc[samp_index,1]
        # testprint
        #print('samplename \t{} \t and corresponding ploidy {} \n'.format(ascatname_list[i], ploidy_file.iloc[samp_index,]))
        sampmatrix = make_sampmatrix(probe_file, ascat_file, ascatname_list[i], matrix_type, tumploidy)
        # list with sample matrices
        fin_list.append(sampmatrix)
        i += 1
        
    # merge into final matrix
    finmatrix = reduce(lambda x, y: pd.merge(x, y, on = ['ProbeID','Chr','Position']), fin_list)
    #tune_finmatrix(finmatrix)
    return finmatrix

#%% test
'''
testlist = ascatfile_list[0:1]
m = make_finmatrix(probe_file,testlist,ascatname_list,ploidy_file,matrix_type="total")
m.info()
#m.loc[m['Chr']=="1"]
m["Chr"].unique()


# the y chromosome probes are probably all NA at the end
'''
##############################################################################
#%% 
''' running code '''

#if __name__ == '__main__': # only execute the code when it is run as a program
    # print out the settings the program is run with
    #print('The program is run with the settings:\ninput genome file:\t{}\noutput file name:\t{}\nminimum ORF length:\t{}\nShine-Dalgarno-seq.:\t{}\nstart codons:\t\t{}\nstop codons:\t\t{}\nalign. error threshold:\t{}\nouput fasta:\t\t{}\n'.format(input_file, output_file, min_length, SD_seq, start_codons, stop_codons, allowed_errors, fasta_output))
total = make_finmatrix(probe_file,ascatfile_list,ascatname_list,ploidy_file,matrix_type="total")
total.to_csv("~/Desktop/MTP_project/Data/SCAN_B/WGS/total_matrix.tsv", sep='\t',encoding='utf-8', index=False)
minor = make_finmatrix(probe_file,ascatfile_list,ascatname_list,ploidy_file,matrix_type="minor")
minor.to_csv("~/Desktop/MTP_project/Data/SCAN_B/WGS/minor_matrix.tsv", sep='\t',encoding='utf-8', index=False)
major = make_finmatrix(probe_file,ascatfile_list,ascatname_list,ploidy_file,matrix_type="major")
major.to_csv("~/Desktop/MTP_project/Data/SCAN_B/WGS/major_matrix.tsv", sep='\t',encoding='utf-8', index=False)
gainloss = make_finmatrix(probe_file,ascatfile_list,ascatname_list,ploidy_file,matrix_type="gainloss")
gainloss.to_csv("~/Desktop/MTP_project/Data/SCAN_B/WGS/gainloss_matrix.tsv", sep='\t',encoding='utf-8', index=False)
amp = make_finmatrix(probe_file,ascatfile_list,ascatname_list,ploidy_file,matrix_type="amplification")
amp.to_csv("~/Desktop/MTP_project/Data/SCAN_B/WGS/amplification_matrix.tsv", sep='\t',encoding='utf-8', index=False)
loh = make_finmatrix(probe_file,ascatfile_list,ascatname_list,ploidy_file,matrix_type="LOH")
loh.to_csv("~/Desktop/MTP_project/Data/SCAN_B/WGS/loh_matrix.tsv", sep='\t',encoding='utf-8', index=False)
#%% the end