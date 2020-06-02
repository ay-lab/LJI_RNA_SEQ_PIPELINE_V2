import pandas as pd
import argparse,os,json,sys
from glob import glob

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-json', '--json_file', type=str)
    opts = parser.parse_args()
    conf_file = opts.json_file
    
    # create all necessary folders
    for folder in [f'2.Internal_files/{folder}' for folder in  ['fastq_input','fastq_filtered','fastp_output','bam_aligned']] + [f'3.QC_files/{folder}' for folder in  ['bed_wiggle','fastp_report','qualimap','qualimap_bamqc']] + [f'4.Output/{folder}' for folder in ['counts','QC_plots','QC_report_history']]:
        os.system(f'mkdir -p {folder}') 

    dict_conf = read_json(conf_file)
    run_file = dict_conf['config']['run_file']
    run_ID = dict_conf['config']['run_ID']
           
    
    # Check if sample_run file exist
    try:
        df_samples = pd.read_csv(run_file,index_col = 0)
        dict_samples = df_samples[run_ID].to_dict()
        
    except KeyError:
        sys.exit('Run does not exist in run sample file!')
    except FileNotFoundError:
        sys.exit('sample run file does not exist!')
        
    # Find samples different from previous run
    try:
        with open('4.Output/lastest_run.txt','r') as f:
            previous_run = f.readline().replace('\n','')
        run_sample_list = list(df_samples[df_samples[previous_run]!=df_samples[run_ID]].index) 

    except FileNotFoundError:
        print('No previous run found, rerun all samples!')
        run_sample_list = list(df_samples.index)
    
    # move fastq files if something changed
    if len(run_sample_list) > 0:
        for sample_ID in run_sample_list:
            fastq_f = []
            fastq_r = []
            for path in dict_samples[sample_ID].split(','):
                fastq_f.append(f'1.Fastq_input/{path}/{sample_ID}_R1.fastq.gz')
                fastq_r.append(f'1.Fastq_input/{path}/{sample_ID}_R2.fastq.gz') 
            os.system(f"cat {' '.join(fastq_f)} > 2.Internal_files/fastq_input/{sample_ID}_R1.fastq.gz")
            os.system(f"cat {' '.join(fastq_r)} > 2.Internal_files/fastq_input/{sample_ID}_R2.fastq.gz")
     
    os.system(f"echo {run_ID} >  4.Output/lastest_run.txt")            
                      
def read_json(file = 'conf_RNA_Seq.json'):
    with open(file) as json_file:
        conf = json.load(json_file)
    return conf                      
                      
if __name__ == "__main__":
    main()