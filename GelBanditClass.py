import os
import subprocess
import pandas as pd
import papermill as pm
from glob import glob
import plotly.express as px
from lxml import etree
from pathlib import Path
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
import shutil
import copy
from Bio import SeqIO
from os.path import basename
import time
import argparse
import sys
import re

class GelBandit:
    def __init__(self, MQ_version, MQ_path, db_path,
                 gel_notebook_template, chromo_notebook_template,
                 fasta_folder, output_folder, conditions, db, reattribute,user_input_params,
                 prot_quantification="Razor + Unique", num_missed_cleavages=2,
                 id_parse_rule=">.*\|(.*)\|",desc_parse_rule=">(.*)",andromeda_path="C:\Temp\Andromeda",
                 fixed_mods=["Carbamidomethyl (C)"], enzymes=["Trypsin/P"], use_enzyme_first_search_str = "True",
                 fs_enzymes= ["Trypsin/P"], var_mods = ["Oxidation (M)", "Acetyl (Protein N-term)"],
                 second_peptide_str = "False", match_between_runs = "False", num_threads = 16,
                 MQ_params=None,
                 extra_POIs=None,
                 intensity_threshold=40, ppm_tolerance=500, rt_match_window=5 
                 ):
        
        self.MQ_path = MQ_path
        self.MQ_version = MQ_version
        self.validate_MQ_path()
        self.database_path = db_path
        self.chromo_notebook_template = chromo_notebook_template
        self.gel_notebook_template = gel_notebook_template
        
        self.fasta_folder = fasta_folder
        self.output_folder = output_folder
        self.conditions = conditions
        self.db = db
        self.reattribute = reattribute
        # Optional args
        # MaxQuant options
        self.prot_quantification = prot_quantification
        self.num_missed_cleavages =  num_missed_cleavages
        self.fixed_mods = fixed_mods
        self.enzymes = enzymes
        self.use_enzyme_first_search_str = use_enzyme_first_search_str
        self.fs_enzymes = fs_enzymes
        self.var_mods = var_mods
        self.second_peptide_str = second_peptide_str
        self.match_between_runs = match_between_runs
        self.num_threads = num_threads
        self.id_parse_rule = id_parse_rule
        self.desc_parse_rule = desc_parse_rule
        self.andromeda_path = andromeda_path
        # Advanced
        self.user_input_params = user_input_params
        self.MQ_params = MQ_params
        # Extra POIs
        self.extra_POIs = extra_POIs
        # MS chromatogram bp identifier display options
        self.intensity_threshold = intensity_threshold
        self.ppm_tolerance = ppm_tolerance
        self.rt_match_window = rt_match_window

        # Hardcoded args
        #self.colours=['#FF5900', '#DA16FF', '#1CA71C', '#2E91E5', '#E15F99', '#25B0B0',
        # '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1',
        # '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1',
        # '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']

        self.colours=['#FF5900', '#FFBA08', '#F48C06', '#9D0208', '#E15F99', '#25B0B0',
         '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1',
         '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1',
         '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']

        # Define file extensions
        self.fasta_extensions = ['.fa', '.FA', '.fasta', '.FASTA', '.faa', '.FAA']
        
        self.results_op = self.make_results_folder()
        self.sentinel_file_path = Path(self.output_folder) / "stop_requested.sentinel"
        self.conditions_dict = self.make_conditions_dict()
        self.raw_folder = self.get_common_parent_directory()
        self.species_dict = self.load_species_dict()
        self.database_path = self.get_species_filepath()
        self.master_fasta = Path(self.fasta_folder) / 'merged_fasta_folder' / f'custom_fasta_merged_w_{self.db.upper()}_database.fasta'
        self.MQ_op_folder = Path(self.raw_folder) / "combined/"
        self.fasta_files = self.get_files_with_extensions()
        self.POIs = self.get_POIs()

        self.msconvert_path = self.get_msconvert_path()

    def validate_MQ_path(self):
        if not os.path.exists(self.MQ_path):
            raise FileNotFoundError(f"The provided MaxQuant path {self.MQ_path} does not exist.")
        if not os.path.isfile(self.MQ_path):
            raise ValueError(f"The provided MaxQuant path {self.MQ_path} is not a file.")

    def make_results_folder(self):
        # Create results folder if it doesn't exist
        results_op = Path(self.output_folder) / "CNIO_prot_core_results"
        if not os.path.exists(results_op):
            os.makedirs(results_op)
        return Path(results_op)

    def sentinel_file_exists(self):
        return os.path.exists(self.sentinel_file_path)

    def make_conditions_dict(self):
        try:
            ext = os.path.splitext(str(Path(self.conditions)))[1]
            if ext == ".xlsx":
                conditions_df = pd.read_excel(Path(self.conditions), sheet_name="Sheet1")
            elif ext == ".txt" or ext == ".tsv":
                conditions_df = pd.read_csv(Path(self.conditions), sep="\t")
        except Exception as e:
            print(f"Error: {e}. Please use a file in excel or tab-delimited format.")
            return None
        
        # Add a 'Basename' column to the dataframe
        print(conditions_df)
        conditions_df['Basename'] = conditions_df['Raw file path'].apply(lambda x: os.path.basename(str(Path(x))))
        
        # Convert the dataframe to a dictionary with 'raw_file' as keys
        conditions_dict = {}
        for _, row in conditions_df.iterrows():
            raw_file = row['Basename']
            conditions_dict[raw_file] = row.drop('Basename').to_dict()
        return conditions_dict

    # Define function to map species to file location
    def load_species_dict(self):
        # Load the file into a dataframe
        df = pd.read_excel(self.database_path, sheet_name="Sheet1")
    
        # Convert species names to lowercase and create dictionary
        species_dict = {row['Species'].upper(): row['Local path'] for _, row in df.iterrows()}
    
        return species_dict

    def get_species_filepath(self):
        # Convert input to uppercase and fetch the filepath
        species = self.db.upper()
        return Path(self.species_dict.get(species, None))

    def get_common_parent_directory(self):
        parent_dirs = set()
        file_paths = [values['Raw file path'] for values in self.conditions_dict.values()]
        for file_path in file_paths:
            # get the parent directory name
            parent_dir = Path(file_path).parent
            parent_dirs.add(parent_dir)
            print(parent_dir)

        # check if all parent directory names are the same
        if len(parent_dirs) == 1:
            return parent_dirs.pop()
        else:
            raise ValueError("Multiple parent directories found!")

    def get_files_with_extensions(self):
        files = []
        folder=Path(self.fasta_folder)
        for extension in self.fasta_extensions:
            files.extend(folder.glob(f'*{extension}'))
        files=list(set(files))
        return files

    def get_POIs(self):
        # Initialize all_POIs with any provided extra POIs or an empty list
        all_POIs = self.extra_POIs if self.extra_POIs else []
        
        for fasta_file in self.fasta_files:
            with open(fasta_file) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    try:
                        POI = record.id.split('|')[1]
                        all_POIs.append(POI)
                    except IndexError as e:
                        print("Look like one of your custom fasta files isn't Uniprot-formatted.")
                        print(f"{e}")
        
        return all_POIs

    def get_msconvert_path(self):
        msconvert_path = None
        for root_folder in ['C:\\Users','C:\\Program Files', 'C:\\Program Files (x86)']:
            for root, dirs, files in os.walk(root_folder):
                if 'msconvert.exe' in files:
                    msconvert_path = os.path.join(root, 'msconvert.exe')
                    break  # Stop the inner loop if we find msconvert.exe
            if msconvert_path:  # Stop the outer loop if we find msconvert.exe
                break
        print(msconvert_path)
        return msconvert_path

    # Non-init methods
    def check_and_terminate_if_sentinel_exists(self):
        if self.sentinel_file_exists():
            try:
                print("Trying cleanup")
                if hasattr(self, 'results_op'):
                    if Path(self.results_op).exists():
                        print("removing folder")
                        shutil.rmtree(Path(self.results_op))
                if hasattr(self, 'temp_MQ_params'):
                    if Path(self.temp_MQ_params).exists():
                        print("removing temp params")
                        Path(self.temp_MQ_params).unlink()
                    if Path(str(self.temp_MQ_params).split('_temp.xml')[0] + '_updated.xml').exists():
                        print("removing updated params")
                        Path(str(self.temp_MQ_params).split('_temp.xml')[0] + '_updated.xml').unlink()

                # Removing any UNFINISHED MQ results
                # Completed runs will remain in the output folder
                if hasattr(self, 'MQ_op_results'):
                    if Path(self.MQ_op_results).exists():
                        run_times_file = Path(self.MQ_op_results) / 'proc' / '#runningTimes.txt'
                        try:
                            with open(run_times_file, 'r') as f:
                                lines = f.readlines()
                                # Get the last line
                                last_line = lines[-1].strip()
                                # Check if the last line contains the string "Finish writing tables"
                                if "Finish writing tables" not in last_line:
                                    shutil.rmtree(Path(self.MQ_op_results))             
                        except FileNotFoundError:
                            if Path(self.MQ_op_results) / 'txt' / 'allPeptides.txt' == False:
                                shutil.rmtree(Path(self.MQ_op_results))
                if hasattr(self, 'paramterized_nb'):
                    if Path(self.paramterized_nb).exists():
                        Path(self.paramterized_nb).unlink()
                    if Path(self.paramterized_nb.replace('.ipynb','.html')).exists():
                        Path(self.paramterized_nb.replace('.ipynb','.html')).unlink()
                # Removing ALL completed file converts for now if this point of the program was reached
                if hasattr(self, 'converted_op'):
                    if Path(self.converted_op).exists():
                        shutil.rmtree(Path(self.converted_op))
                if hasattr(self, 'paramterized_chromo_nb'):
                    if Path(self.paramterized_chromo_nb).exists():
                        Path(self.paramterized_chromo_nb).unlink()
                    if Path(self.paramterized_chromo_nb.replace('.ipynb','.html')).exists():
                        Path(self.paramterized_chromo_nb.replace('.ipynb','.html')).unlink()
                # Clean up any resources if necessary
                # e.g. close open files, terminate child processes, etc.
                # Checking all files produced in script
                print("Sentinel file detected. Terminating script.")
                exit()  # This will terminate the current script
            except Exception as e:
                print(f"An error occurred during cleanup: {e}")
            finally:
                print("Sentinel file detected. Terminating script.")
                exit()  # This will terminate the current script
    
    def concatenate_fasta_files(self):
        def append_file(inp_file, output_file):
            file_wo_ext,ext = os.path.splitext(inp_file)
            copied_file=Path(f'{file_wo_ext}_temp.fasta')
            shutil.copyfile(inp_file, copied_file)
            with open(copied_file, "a") as file:
                file.write('\n')
                file.close()
            subprocess.run([r'C:\Users\lwoods\AppData\Local\mambaforge\envs\prot_plots_env\Library\usr\bin\cat.exe',str(copied_file),'>>',str(output_file)],shell=True)
            os.remove(copied_file)   
        Path(self.master_fasta).parent.mkdir(parents=True, exist_ok=True)
        if not os.path.exists(self.master_fasta):
            for fasta_file in self.fasta_files:
                append_file(fasta_file, self.master_fasta)
            append_file(self.db, self.master_fasta)
            with open(self.master_fasta, 'r') as infile, open(Path(str(self.master_fasta).replace('.fasta', '_no_empty_lines.fasta')), 'w') as outfile:
                for line in infile:
                    if line.strip():
                        outfile.write(line)
            os.remove(Path(self.master_fasta))
            shutil.move(Path(str(self.master_fasta).replace('.fasta','_no_empty_lines.fasta')), Path(self.master_fasta))   
        else:
            print("Merged fasta already exists--skipping\n")

    def create_MaxQuant_par(self):
        # Create an empty params file for editing if one isn't provided
        if self.MQ_params:
            pass
        else:
            print(f"Creating an MQ par file!")
            new_MQ_params_name = f"mqpar_version_{self.MQ_version}.xml"
            if (Path(self.output_folder) / new_MQ_params_name).exists():
                print("Deleting existing default mqpar")
                (Path(self.output_folder) / new_MQ_params_name).unlink()
            
            self.MQ_params = Path(self.output_folder) / new_MQ_params_name
            process_0 = subprocess.Popen(
            [self.MQ_path, '--create',  str(self.MQ_params)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            while process_0.poll() is None:
                self.check_and_terminate_if_sentinel_exists()
                time.sleep(2)

            result_stdout, result_stderr = process_0.communicate()

            if process_0.returncode != 0:
                print("MaxQuant --create encountered an error.")
                print("Error Output:")
                print(result_stderr)
                return  # Exit the method if the first subprocess fails

            print("MaxQuant --create ran successfully.")
            print("Standard Output:")
            print(result_stdout)

    def edit_MQ_par(self):
        print(Path(self.MQ_params))
        tree = etree.parse(Path(self.MQ_params))
        root = tree.getroot()

        # Editing fasta location
        fastaFile_block = root.find("./fastaFiles")

        for child in list(fastaFile_block)[1:]:
            fastaFile_block.remove(child)

        FastaFileInfo = fastaFile_block.find("FastaFileInfo")
        fasta_path = FastaFileInfo.find('fastaFilePath')
        fasta_path.text = os.path.basename(self.master_fasta)

        # Editing raw files and corresponding attributes
        branch_names = ["./filePaths", "./experiments", "./fractions", "./ptms", "./paramGroupIndices", "./referenceChannel"]
        branch_list = [root.find(branch_name) for branch_name in branch_names]

        for branch in branch_list:
            for child in list(branch)[1:]:
                branch.remove(child)

        child_tags = ['string', 'string', 'short', 'boolean', 'int', 'string']
        child_list = [branch.find(tag) for tag, branch in zip(child_tags, branch_list)]

        ET_raw_files = list(self.conditions_dict.keys())
        experiments_list = [details['Experiment'] for details in self.conditions_dict.values()]
        common_length = len(self.conditions_dict)
        fractions = ['32767'] * common_length
        ptms = ['False'] * common_length
        pGIs = ['0'] * common_length
        refChannels = [''] * common_length

        par_text_lists = [ET_raw_files, experiments_list, fractions, ptms, pGIs, refChannels]

        for i, new_par_texts in enumerate(zip(*par_text_lists)):
            for par_child, new_par_text in zip(child_list, new_par_texts):
                if i > 0:
                    new_dupe = copy.deepcopy(par_child)
                    new_dupe.text = new_par_text
                    branch_list[child_list.index(par_child)].append(new_dupe)
                else:
                    par_child.text = new_par_text

        # Only if not using input MQ params
        if self.user_input_params == False:
            # Quantification type
            quant_mode_branch = root.find("./quantMode")
            quant_mode_dict = {'All': '0', 'Razor + Unique': '1', 'Unique': '2'}
            quant_mode_branch.text = quant_mode_dict[self.prot_quantification]

            # Edits to all parameter groups
            for param_group in root.findall("./parameterGroups/parameterGroup"):

                # Max missed cleavages
                missed_cleavage_branch = param_group.find("maxMissedCleavages")
                if missed_cleavage_branch is not None:
                    missed_cleavage_branch.text = str(self.num_missed_cleavages)

                # Fixed modifications
                ## self.fixed_mods
                fixed_mod_branches = param_group.findall("fixedModifications")
                for fixed_mod_branch in fixed_mod_branches:
                    # Remove existing fixed modifications if needed
                    for child in list(fixed_mod_branch):
                        fixed_mod_branch.remove(child)
                    # Add new fixed modifications
                    for new_fixed_mod in self.fixed_mods:
                        new_fixed_mod_elem = etree.SubElement(fixed_mod_branch, "string")
                        new_fixed_mod_elem.text = new_fixed_mod

                # Enzymes
                ## self.enzymes(list)
                enzymes_branches = param_group.findall("enzymes")
                for enzymes_branch in enzymes_branches:
                    # Remove existing fixed modifications if needed
                    for child in list(enzymes_branch):
                        enzymes_branch.remove(child)
                    # Add new fixed modifications
                    for new_enzyme in self.enzymes:
                        new_enzyme_elem = etree.SubElement(enzymes_branch, "string")
                        new_enzyme_elem.text = new_enzyme

                # Use enzymes first search
                ## self.use_enzyme_first_search_str
                use_enzyme_first_search_branch = param_group.find("useEnzymeFirstSearch")
                if use_enzyme_first_search_branch is not None:
                    use_enzyme_first_search_branch.text = str(self.use_enzyme_first_search_str)

                # First search enzymes 
                ## self.fs_enzymes(list)
                fs_enzymes_branches = param_group.findall("enzymesFirstSearch")
                for fs_enzymes_branch in fs_enzymes_branches:
                    # Remove existing fixed modifications if needed
                    for child in list(fs_enzymes_branch):
                        fs_enzymes_branch.remove(child)
                    # Add new fixed modifications
                    for new_fs_enzyme in self.fs_enzymes:
                        new_fs_enzyme_elem = etree.SubElement(fs_enzymes_branch, "string")
                        new_fs_enzyme_elem.text = new_fs_enzyme

                # Var modifications
                ## self.var_mods
                var_mod_branches = param_group.findall("variableModifications")
                for var_mod_branch in var_mod_branches:
                    # Remove existing fixed modifications if needed
                    for child in list(var_mod_branch):
                        var_mod_branch.remove(child)
                    # Add new fixed modifications
                    for new_var_mod in self.var_mods:
                        new_var_mod_elem = etree.SubElement(var_mod_branch, "string")
                        new_var_mod_elem.text = new_var_mod


            
            # Editing parse rules
            identifierParseRule = FastaFileInfo.find("identifierParseRule")
            identifierParseRule.text = str(self.id_parse_rule)
            descriptionParseRule = FastaFileInfo.find("descriptionParseRule")
            descriptionParseRule.text = str(self.desc_parse_rule)

            #Editing fixed search folder for Andromeda
            fixedSearchFolder = root.find("./fixedSearchFolder")
            fixedSearchFolder.text = str(self.andromeda_path)

            # Editing option to try second peptide
            ## self.second_peptide_str
            secondPeptide = root.find("./secondPeptide")
            secondPeptide.text = str(self.second_peptide_str)

            # Editing option to match between runs
            ## self.match_between_runs
            matchBetweenRuns = root.find("./matchBetweenRuns")
            matchBetweenRuns.text = str(self.match_between_runs)

            ## self.num_threads
            numThreads = root.find("./numThreads")
            numThreads.text = str(self.num_threads)

            

        et = etree.ElementTree(root)
        self.temp_MQ_params = Path(self.output_folder, os.path.basename(str(Path(self.MQ_params))).split('.xml')[0] + '_temp.xml')
        et.write(self.temp_MQ_params, pretty_print=True)

    def run_MaxQuant(self):
        win_MQ_params_updated = Path(str(self.temp_MQ_params).split('_temp.xml')[0] + '_updated.xml')
        print(win_MQ_params_updated)
        master_folder = str(Path(self.master_fasta).parent)
        print(master_folder)
        print(self.MQ_op_folder)

        if not self.MQ_op_folder.exists():

            self.check_and_terminate_if_sentinel_exists()
            
            if win_MQ_params_updated.exists():
                win_MQ_params_updated.unlink()

            # Start changeFolder subprocess
            process_1 = subprocess.Popen(
                [self.MQ_path, str(self.temp_MQ_params), '--changeFolder', str(win_MQ_params_updated), str(master_folder), str(self.raw_folder)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Periodically check if the sentinel file exists while the subprocess is running
            while process_1.poll() is None:
                self.check_and_terminate_if_sentinel_exists()
                time.sleep(2)  # You can adjust this sleep duration

            # Fetch the results after the process ends
            result_stdout, result_stderr = process_1.communicate()

            if process_1.returncode == 0:
                print("MaxQuant changeFolder ran successfully.")
                print("Standard Output:")
                print(result_stdout)
            else:
                print("MaxQuant changeFolder encountered an error.")
                print("Error Output:")
                print(result_stderr)
            
            # Start search subprocess
            process_2 = subprocess.Popen(
                [self.MQ_path, str(win_MQ_params_updated)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Again, periodically check if the sentinel file exists while the subprocess is running
            while process_2.poll() is None:
                self.check_and_terminate_if_sentinel_exists()
                time.sleep(2)

            # Fetch the results after the process ends
            MQ_result_stdout, MQ_result_stderr = process_2.communicate()

            if process_2.returncode == 0:
                print("MaxQuant search ran successfully.")
                print("Standard Output:")
                print(MQ_result_stdout)
            else:
                print("MaxQuant search encountered an error.")
                print("Error Output:")
                print(MQ_result_stderr)
        else:
            print("MaxQuant output folder already available.")

    def get_poi_colour_map(self):
        colours_in_use = self.colours[0:len(self.POIs)]
        poi_colours_map = dict(zip(sorted(self.POIs), colours_in_use))
        return poi_colours_map
   
    def generate_excel_results(self):
        self.check_and_terminate_if_sentinel_exists()
        self.poi_colours_map = self.get_poi_colour_map()
        print('Generating excel files')
        print(self.POIs)
        prot_group_file = Path(self.MQ_op_folder / 'txt/proteinGroups.txt')
        peptides_file = Path(self.MQ_op_folder / 'txt/peptides.txt')
        msms_file = Path(self.MQ_op_folder / 'txt/msms.txt')
        prot_groups_df = pd.read_csv(prot_group_file, sep='\t',low_memory=False)
        peptides_df=pd.read_csv(peptides_file, sep='\t',low_memory=False)
        msms_df=pd.read_csv(msms_file, sep='\t',low_memory=False)

        base_cols = ['Majority protein IDs', 'Protein names', 'Gene names', 'Fasta headers', 'Intensity']
        experiments = [values['Experiment'] for values in self.conditions_dict.values()]
        experiment_cols = [f'Peptides {exp}' for exp in experiments] + [f'Sequence coverage {exp} [%]' for exp in experiments] + [f'Intensity {exp}' for exp in experiments]
        summary_prot_group_cols = base_cols + experiment_cols + ['Potential contaminant', 'Score']

        def extract_gene_names(header):
            if ";" in str(header):
                headers = str(header).split(";")
            else:
                headers = [str(header)]
            gene_names = [re.search(r'GN=([\w\d]+)', h).group(1) if re.search(r'GN=([\w\d]+)', h) else '' for h in headers]
            return ";".join(gene_names)

        def extract_protein_names(header):
            if ";" in str(header):
                headers = str(header).split(";")
            else:
                headers = [str(header)]
            protein_names = [re.search(r'\|([^|]+)\s+OS=', h).group(1) if re.search(r'\|([^|]+)\s+OS=', h) else '' for h in headers]
            return ";".join(protein_names)
        
        if 'Gene names' not in prot_groups_df.columns:
            prot_groups_df['Gene names'] = prot_groups_df['Fasta headers'].apply(extract_gene_names)

        if 'Protein names' not in prot_groups_df.columns:
            prot_groups_df['Protein names'] = prot_groups_df['Fasta headers'].apply(extract_protein_names)

        # Perform reattribution, if required
        # If re-attributing, colours will be per POI
        # Else, colours will be per Majority protein ID
        def reattribute_proteingroups(df, pois):
            df['Reattributed Protein'] = df['Majority protein IDs']

            for idx, row in df.iterrows():
                proteins_str = row['Majority protein IDs']
                if not isinstance(proteins_str, str):
                    continue  # Skip rows without a valid string value
                proteins = set(proteins_str.split(';'))
                overlap = proteins.intersection(pois)
                
                if not overlap:
                    continue
                
                if row['Majority protein IDs'] in overlap:
                    df.at[idx, 'Reattributed Protein'] = row['Majority protein IDs']
                else:
                    df.at[idx, 'Reattributed Protein'] = ';'.join(overlap)

            return df

        def reattribute_peptides(df, pois):
            df['Reattributed Protein'] = df['Proteins']
            for idx, row in df.iterrows():
                proteins_str = row['Proteins']
                if not isinstance(proteins_str, str):
                    continue  # Skip rows without a valid string value
                proteins = set(proteins_str.split(';'))
                overlap = proteins.intersection(pois)
                if not overlap:
                    continue
                leading_majority = row.get('Leading razor protein', None)
                if leading_majority and leading_majority in overlap:
                    df.at[idx, 'Reattributed Protein'] = leading_majority
                else:
                    df.at[idx, 'Reattributed Protein'] = ';'.join(overlap)
            return df
        
        if self.reattribute==True:
            prot_groups_df=reattribute_proteingroups(prot_groups_df, self.POIs)
            peptides_df=reattribute_peptides(peptides_df, self.POIs)
            base_cols.append('Reattributed Protein')
            prot_cols_list=['Reattributed Protein']
            pept_cols_list=['Reattributed Protein']
            summary_prot_group_cols.append('Reattributed Protein')
        else:
            prot_cols_list=['Majority protein IDs','Proteins IDs']
            pept_cols_list=['Leading razor protein','Proteins']
            
        # Sort the DataFrames
        prot_groups_df = prot_groups_df.sort_values(by='Intensity', ascending=False)
        summary_prot_df= prot_groups_df[summary_prot_group_cols]

        def get_row_color(row, cols_list, poi_colours_map):
            proteins = set()
            for col_name in cols_list:
                try:
                    if pd.notna(row[col_name]):
                        proteins.update(row[col_name].split(';'))
                except KeyError:
                    continue  # Ignore columns that don't exist in the row

            # Priority-based coloring
            for poi, colour in poi_colours_map.items():
                if poi in proteins:
                    return colour

            return None  # Default no-color
        
        # Define row highlight function
        def calculate_brightness(hex_color):
            # Calculate brightness based on RGB values
            r, g, b = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)
            brightness = (r * 299 + g * 587 + b * 114) / 1000
            return brightness

        def get_text_color(hex_color):
            # Determine whether to use white or black text based on brightness
            brightness = calculate_brightness(hex_color)
            return "#000000" if brightness > 128 else "#FFFFFF"
        
        def highlight_rows(row, cols_list, poi_colours_map):
            colour = get_row_color(row, cols_list, self.poi_colours_map)
            if colour:
                text_color = get_text_color(colour)
                style = f'background-color: {colour}; color: {text_color}'
                return [style] * len(row)
            return [''] * len(row)

        self.excel_op = Path(self.results_op / "MaxQuant_protein_groups_results.xlsx")
        # Write the DataFrames to an Excel file with style and sheets
        with pd.ExcelWriter(self.excel_op) as writer:
            self.check_and_terminate_if_sentinel_exists()
            #prot_groups_df.style.apply(highlight_rows, args=('Majority protein IDs,Proteins IDs', POIs, colours), axis=1).to_excel(writer, sheet_name='Full_proteinGroups', index=False)
            prot_groups_df.style.apply(highlight_rows, args=(prot_cols_list, self.poi_colours_map), axis=1).to_excel(writer, sheet_name='Full_proteinGroups', index=False)
            summary_prot_df.style.apply(highlight_rows, args=(prot_cols_list, self.poi_colours_map), axis=1).to_excel(writer, sheet_name='Summary_proteinGroups', index=False)
            peptides_df.style.apply(highlight_rows, args=(pept_cols_list, self.poi_colours_map), axis=1).to_excel(writer, sheet_name='Peptides', index=False)
            msms_df.style.apply(highlight_rows, args=(['Proteins'], self.poi_colours_map), axis=1).to_excel(writer, sheet_name='msms', index=False)

    def jupyter_conversion(self, parameterized_nb):
        # Must be in the right environment!
        output_html = Path(str(parameterized_nb).replace('ipynb', 'html'))
        print(output_html)
        subprocess.run(["jupyter", "trust", str(parameterized_nb)], capture_output=True, text=True)
        subprocess.run(["jupyter", "nbconvert", "--no-input", "--to", "html", str(parameterized_nb), str(output_html)], capture_output=True, text=True)

    def run_parameterized_nb(self):
        self.check_and_terminate_if_sentinel_exists()
        self.op_notebook = Path(self.output_folder) / "Gel_bandit_plotter.ipynb"
        print(self.gel_notebook_template)
        print(self.op_notebook)
        print("Arguments:")
        print(str(Path(self.excel_op)))
        print(str(self.master_fasta))
        print(self.POIs)
        experiments = [values['Experiment'] for values in self.conditions_dict.values()]
        print(experiments)
        print(self.prot_quantification)
        pm_exec = pm.execute_notebook(
            str(self.gel_notebook_template),
            output_path=str(self.op_notebook),
            parameters=dict(
                excel_path=Path(self.excel_op).as_posix(),
                fasta_file=Path(self.master_fasta).as_posix(),
                POIs=self.POIs,
                experiments=experiments,
                identification=self.prot_quantification,
                colours=self.colours,
                reattribute_peptides=self.reattribute
                )
        )
        self.jupyter_conversion(self.op_notebook)

    def convert_files(self):
        if self.msconvert_path:
            
            self.check_and_terminate_if_sentinel_exists()
            
            self.converted_op = Path(self.output_folder) / "mzML_files"
            if not os.path.exists(self.converted_op):
                os.makedirs(self.converted_op)

            raw_files = [details['Raw file path'] for details in self.conditions_dict.values()]

            if all([os.path.exists(raw_file) for raw_file in raw_files]):
                print("All raw files exist.")
            else:
                print("Not all raw files exist!")
                return
            
            # Determine the expected mzML filenames
            self.mzML_files = [str(self.converted_op / str(os.path.basename(os.path.splitext(raw_file)[0])+ '.mzML')) for raw_file in raw_files]

            # Check if mzML files already exist
            if all([os.path.exists(mzML_file) for mzML_file in self.mzML_files]):
                print("All mzML files already exist. No conversion performed.")
                return

            else:

                # If they don't all exist, proceed with the conversion
                BASE_COMMAND = [
                    "--zlib",
                    "--outdir", f'{str(self.converted_op)}',
                    "--filter", "peakPicking vendor msLevel=1-",
                    "--filter", "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState> File:\"\"\"^<SourcePath^>\"\"\", NativeID:\"\"\"^<Id^>\"\"\""
                ]

                raw_files = [str(Path(raw_file)) for raw_file in raw_files]
                command = [str(Path(self.msconvert_path))] + BASE_COMMAND + raw_files

                # Start search subprocess
                process_3 = subprocess.Popen(
                    command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                    )

                # Again, periodically check if the sentinel file exists while the subprocess is running
                while process_3.poll() is None:
                    self.check_and_terminate_if_sentinel_exists()
                    time.sleep(2)

                # Fetch the results after the process ends
                MSConvert_result_stdout, MSConvert_result_stderr = process_3.communicate()

                if process_3.returncode == 0:
                    print("MSConvert ran successfully.")
                    print("Standard Output:")
                    print(MSConvert_result_stdout)
                else:
                    print("MSConvert encountered an error.")
                    print("Error Output:")
                    print(MSConvert_result_stderr)        

    def run_parameterized_chromo_nb(self):
        self.check_and_terminate_if_sentinel_exists()
        self.op_chromo_notebook = Path(self.output_folder) / f"Chromato_plotter_parameterized.ipynb"
        print(self.chromo_notebook_template)
        print(self.op_chromo_notebook)
        print("Arguments:")
        pm_exec = pm.execute_notebook(
            str(self.chromo_notebook_template),
            output_path=str(self.op_chromo_notebook),
            parameters=dict(
                MQ_op_folder=str(Path(self.MQ_op_folder)),
                mzml_files=self.mzML_files,
                intensity_thresh = self.intensity_threshold,
                ppm_tolerance = self.ppm_tolerance,
                rt_match_window = self.rt_match_window,
                POIs = self.POIs
            )
        )
        self.jupyter_conversion(self.op_chromo_notebook)


def comma_separated_string_to_list(comma_separated_string):
    return comma_separated_string.split(',')

def str_to_bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == "__main__":

# FOR TESTING! ###
    sys.argv = [
        "GelBanditClass",  # normally this is the script name, it's not important here
        "--MQ_version", "2.1.4.0",
        "--MQ_path", r"C:\Users\lwoods\Desktop\MaxQuant_v2.1.4.0 (1)\MaxQuant v2.1.4.0\bin\MaxQuantCmd.exe",
        "--db_path", r"C:\Users\lwoods\Documents\LW_Projects_folder\general\database_map.xlsx",
        "--gel_bandit_plotter", r"C:\Users\lwoods\Documents\Python Scripts\gel_bandit_plotter.ipynb",
        "--chromo_plotter", r"C:\Users\lwoods\Documents\Python Scripts\Chromato-plotter_v4.ipynb",
        "--fasta_folder", r"C:\Users\lwoods\Documents\LW_Projects_folder\Poster_data\GelBandit\fasta_folder",
        "--output_folder", r"C:\Users\lwoods\Documents\LW_Projects_folder\Poster_data\GelBandit",
        "--conditions", r"C:\Users\lwoods\Documents\LW_Projects_folder\Poster_data\GelBandit\conditions.xlsx",
        "--db", "Insect (SF9)",
        "--reattribute", "no",
        "--user_input_params", "False",
        #"--prot_quantification", "Razor + Unique",
        #"--num_missed_cleavages", "2",
        #"--rt_match_window", "2",
        "--extra_POIs", "POI_1"
    ]
##################
   
    parser = argparse.ArgumentParser(description='GelBandit arguments')
    parser.add_argument('--MQ_version', type=str, required=True, help='MaxQuant version.')
    parser.add_argument('--MQ_path', type=str, required=True, help='Path to MaxQuant.')
    parser.add_argument('--db_path', type=str, required=True, help='Database path.')
    parser.add_argument('--gel_bandit_plotter', type=str, required=True, help='GelBandit plotter ipynb template path.')
    parser.add_argument('--chromo_plotter', type=str, required=True, help='Chromatography bp id plotter ipynb template path.')
    parser.add_argument('--fasta_folder', type=str, required=True, help='Path to the fasta folder.')
    parser.add_argument('--output_folder', type=str, required=True, help='Path to the output folder.')
    parser.add_argument('--conditions', type=str, required=True, help='Path to the conditions file.')
    parser.add_argument('--db', type=str, required=True, help='Database choice.')
    parser.add_argument('--reattribute', type=str, required=True, help='Reattribute peptides option (yes or no).')
    parser.add_argument('--user_input_params', type=str_to_bool, required=True, help='User has supplied an MQ params file.')

    # Optional arguments
    parser.add_argument('--prot_quantification', required=False, type=str, default="Razor + Unique", help='Protein quantification method (Default: Razor + Unique).')
    parser.add_argument('--num_missed_cleavages', required=False, type=int, default=2, help='Number of missed cleavages permitted (Default: 2).')

    parser.add_argument('--fixed_mods', required=False, type=list, default=["Carbamidomethyl (C)"], help='default=["Carbamidomethyl (C)"]')
    parser.add_argument('--enzymes', required=False, type=list, default=["Trypsin/P"], help='default=["Trypsin/P"]')
    parser.add_argument('--use_enzyme_first_search_str', required=False, type=str, default="True", help='default=True')
    parser.add_argument('--fs_enzymes', required=False, type=list, default=["Trypsin/P"], help='default=["Trypsin/P"]')
    parser.add_argument('--var_mods', required=False, type=list, default=["Oxidation (M)", "Acetyl (Protein N-term)"], help='default=["Oxidation (M)", "Acetyl (Protein N-term)"]')
    parser.add_argument('--second_peptide_str', required=False, type=str, default="False", help='default=False')
    parser.add_argument('--match_between_runs', required=False, type=str, default="False", help='default=False')
    parser.add_argument('--num_threads', required=False, type=int, default=16, help='Number of threads to use (Default: 16).')
    
    parser.add_argument('--id_parse_rule', required=False, type=str, default=">.*\|(.*)\|", help='Parsing rules for fasta IDs. Changing not recommended! (Default: >.*\|(.*)\| (Uniprot)).')
    parser.add_argument('--desc_parse_rule', required=False, type=str, default=">(.*)", help='Parsing rules for fasta descriptions. Changing not recommended! (Default: >(.*) (Uniprot)).')
    parser.add_argument('--andromeda_path', required=False, type=str, default="C:\Temp\Andromeda", help='Location for storing Andromeda fixed search folder (Default: C:\Temp\Andromeda).')
    parser.add_argument('--MQ_params', required=False, type=str, default=None, help='Path to the MQ params file (Default: None.')
    parser.add_argument('--extra_POIs', required=False, type=comma_separated_string_to_list, default=None, help='Comma-separated extra proteins of interest.')
    parser.add_argument('--intensity_threshold', required=False, type=float, default=40, help='Threshold (%) to display base peak labels (Default: 40).')
    parser.add_argument('--ppm_tolerance', required=False, type=float, default=500, help='Match tolerance for matching base peak m/z values to MaxQuant (Default: 500).')
    parser.add_argument('--rt_match_window', required=False, type=float, default=5, help='Match window (minutes) for matching base peaks to MaxQuant (Default: 5).')
    
    args = parser.parse_args()
    print(args)
    processor = GelBandit(args.MQ_version, args.MQ_path, args.db_path,
                        args.gel_bandit_plotter, args.chromo_plotter, 
                        args.fasta_folder, args.output_folder,
                        args.conditions, args.db, args.reattribute,args.user_input_params,
                        prot_quantification=args.prot_quantification, num_missed_cleavages=args.num_missed_cleavages,
                        id_parse_rule=args.id_parse_rule, desc_parse_rule=args.desc_parse_rule, andromeda_path=args.andromeda_path,
                        fixed_mods=args.fixed_mods, enzymes=args.enzymes, use_enzyme_first_search_str = args.use_enzyme_first_search_str,
                        fs_enzymes=args.fs_enzymes, var_mods=args.var_mods,
                        second_peptide_str=args.second_peptide_str, match_between_runs=args.match_between_runs, num_threads=args.num_threads,
                        MQ_params=args.MQ_params,
                        extra_POIs=args.extra_POIs,
                        intensity_threshold=args.intensity_threshold, ppm_tolerance=args.ppm_tolerance, rt_match_window=args.rt_match_window  
                        )
    processor.concatenate_fasta_files()
    processor.create_MaxQuant_par()
    processor.edit_MQ_par()
    processor.run_MaxQuant()
    processor.generate_excel_results()
    processor.run_parameterized_nb()
    processor.convert_files()
    processor.run_parameterized_chromo_nb()
