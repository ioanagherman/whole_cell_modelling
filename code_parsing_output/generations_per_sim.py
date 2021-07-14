# This file is used to read from a given output folder and write to a df the gene index corresponding to the
# knockdowns that caused cell death. It also writes the generation where the death occurs (the data is not being)
# outputted
# You will have to change the paths to your own directories

import os
import pandas as pd
import datetime
import json

def search_error_in_block(output_path, type_sim):
    # This method will look inside the blocks and launchers directories that fireworks created, find the corresponding ones for each simulation that failed
    # and output what is stored in error file
    block_dirs = os.listdir('/newhome/ig13470/wholecell3/wcEcoli/wholecell/fireworks')
    if type_sim=='kd/ko':
        date = [output_path.split('/')[-5].split('.')[0]]
    elif type_sim=='w':
        date = [output_path.split('/')[7].split('.')[0]]

    date_d = datetime.datetime.strptime(date[0], '%Y%m%d')
    date.append((date_d + datetime.timedelta(days=1)).strftime("%Y%m%d"))
    date.append((date_d - datetime.timedelta(days=1)).strftime("%Y%m%d"))
    date.append((date_d + datetime.timedelta(days=2)).strftime("%Y%m%d"))
    date.append((date_d - datetime.timedelta(days=2)).strftime("%Y%m%d"))

    for block in block_dirs:
        if block[:5]=='block':
            for launcher in os.listdir('/newhome/ig13470/wholecell3/wcEcoli/wholecell/fireworks/' + block):
                date_launch = launcher.split('_')[1].replace('-', '')[:8]
                if date_launch in date:
                    try:
                        with open('/newhome/ig13470/wholecell3/wcEcoli/wholecell/fireworks/' + block + '/' + launcher + "/FW.json", "r") as read_file:
                            d = str(json.load(read_file))
                    except:
                        print('/newhome/ig13470/wholecell3/wcEcoli/wholecell/fireworks/' + block + '/' + launcher + "/FW.json")
                        continue

                    try:
                        if output_path in d:
                            with open('/newhome/ig13470/wholecell3/wcEcoli/wholecell/fireworks/' + block + '/' + launcher + "/FW_job.error") as f:
                                lines = f.readlines()
                                return lines, '/newhome/ig13470/wholecell3/wcEcoli/wholecell/fireworks/' + block + '/' + launcher + "/FW_job.error"
                    except:
                        continue
    return 0, 0


def analyse_generations(path, df, id, type_sim):
    generation_p = []
    list_p = []
    generation_s = []
    list_s = []

    for generations in os.listdir(path):
        if generations.split('_')[0] == 'generation':
            for gen_d in os.listdir(path + '/' + generations):
                if len(os.listdir(path + '/' + generations + '/' + gen_d + '/plotOut')) == 0:
                    list_p.append(id)
                    generation_p.append(generations + '/' + gen_d)

                if len(os.listdir(path + '/' + generations + '/' + gen_d + '/simOut')) == 0:
                    list_s.append(id)
                    generation_s.append(generations + '/' + gen_d)

    if len(list_p) > 0 and len(list_s) > 0 and list_s[0] == list_p[0]:
        error, err_path = search_error_in_block(path + '/' + min(generation_p), type_sim)
        row = [list_p[0].lstrip("0"), path, error, 'N', 'N', err_path, min(generation_s),
               min(generation_p)]
        df = df.append(pd.Series(row, index=list(df.columns)), ignore_index=True)
    elif len(list_p) > 0:
        error, err_path = search_error_in_block(path + '/' + min(generation_p), type_sim)
        row = [list_p[0].lstrip("0"), path, error, 'Y', 'N', err_path, None,
               min(generation_p)]
        df = df.append(pd.Series(row, index=list(df.columns)), ignore_index=True)
    return df

def missing_outputs(out_dir):
    df = pd.DataFrame(columns = ['GeneIdx', 'cell_path', 'Error', 'SimData', 'Plot', 'ErrorPath', 'Generation_death_sim', 'Generation_death_plot'])
    for dirs in os.listdir(out_dir):
        if os.path.isdir(out_dir+dirs):
            for dirs1 in os.listdir(out_dir+dirs):
                if dirs1.split('_')[0][:4]=='gene':
                    # If you are running several initial seeds, you might have to change the '000000' since you will have extra seeds
                    path = out_dir+dirs+'/'+dirs1+'/'+'000000'
                    df = analyse_generations(path, df, dirs1.split('_')[1], 'kd/ko')
    df.to_csv(out_dir + '/sim_with_errors_2.csv', index=False)
    return df


def missing_outputs_wildtype(output_path_wildtype):
    df = pd.DataFrame(columns=['GeneIdx', 'cell_path', 'Error', 'SimData', 'Plot', 'ErrorPath', 'Generation_death_sim', 'Generation_death_plot'])
    for dirs in os.listdir(output_path_wildtype):
        if dirs[:4]=='wild':
            for mcell in os.listdir(output_path_wildtype+'/'+dirs):
                path = output_path_wildtype+'/'+dirs + '/' + mcell + '/'
                df = analyse_generations(path, df, dirs + '/' + mcell, 'w')
    df.to_csv(output_path_wildtype + '/sim_with_errors_2.csv', index=False)
    return df


if __name__ == '__main__':
    #This is the directory where the simulations are stored
    output_path = '/newhome/ig13470/wholecell3/wcEcoli/out/kd_factor_1/'
    print(missing_outputs(output_path))
    #output_path_wildtype = '/newhome/ig13470/wholecell3/wcEcoli/out/wildtype/20210409.233442__WildType_Run_4gen50seed'
    #print(missing_outputs_wildtype(output_path_wildtype))