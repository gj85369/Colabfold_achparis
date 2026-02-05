from alphafold.data.tools import hhsearch
from alphafold.data.tools import hmmsearch
from alphafold.data.tools import jackhmmer
from alphafold.data.tools import hhblits

from alphafold.data import pipeline
from alphafold.data import parsers
import os

import subprocess


class msa_create:
    def __init__(self, name):

        self.name = name
        self.fasta = []    
        self.cwd = os.getcwd()+name
        self.db_list = []  
        self.db_root_path = ''
        self.alphafold_path = ''
        self.jackhmmer_path = ''
        self.hhblits_binary_path='/usr/bin/jackhmmer'   
        self.hhmer_binary_path='/usr/bin/jackhmmer' 
        self.msa_dict = {}
        self.seq_dict = None
        self.mod_db = False
        self.output_path = './'
        self.printout = False
        self.easl_path = None
        self.max_num = None


    def init_db_list(self):
        if self.mod_db == True:
            self.db_list = ['DBDB_mod', 'pdb_seqres']
        else:
            self.db_list = ['uniprot', 'pdb_seqres']


    def init_db_configs(self):
        if self.max_num == None:
            dbmax = 100000
            upmax = 10000
            pmax = 1001
        else:
            dbmax = self.max_num
            upmax = self.max_num
            pmax = self.max_num            
        raw_db_dict = {
            'DBDB_mod': {
                'db_name': 'DBDB_mod',
                'runner': 'hmmer',
                'db_path': self.db_root_path +'/DBDB_mod',
                'num_streamed_chunks': None,
                'z_value': None,
                'max_hits': dbmax,
                'pairable': False},
            'uniprot': {
                'db_name': 'uniprot',
                'runner': 'hmmer',
                'db_path': self.db_root_path +'/uniprot_sprot.fasta',
                'num_streamed_chunks': None,
                'z_value': None,
                'max_hits': upmax,
                'pairable': False},
            'pdb_seqres': {
                'db_name': 'pdbrcsb',
                'runner': 'hmmer',
                'db_path': self.db_root_path + '/pdb_seqres.txt',
                'num_streamed_chunks': None,
                'z_value': None,
                'max_hits': pmax,
                'pairable': False},}
        db_configs = []
        for db in self.db_list:
            db_configs.append(raw_db_dict[db])
        self.db_configs = db_configs
    def make_domtab(self, output_dir, afasta):
        cmd = [self.easl_path, output_dir+'/pdbrcsb.sto', '>' , output_dir+'/domtbl']

        print(cmd)
        ok = subprocess.call(' '.join(cmd), shell=True)
        #ok = subprocess.check_output(cmd)
    def get_MSAs_for_seq_std(self, seq_str,
                             output_dir):
        fasta_path = os.path.join(output_dir, 'query.fasta')
        with open(fasta_path, 'wt') as f:
            f.write(f'>query\n{seq_str}')
        #if self.make_templates:
        #    self.make_domtab(output_dir, fasta_path)
        runner_configs = {
            'hmmer': {
                'runner': jackhmmer.Jackhmmer,
                'binary': self.hhmer_binary_path,
                'format': 'sto',
                'parser': parsers.parse_stockholm},
            'hhblits': {
                'runner': hhblits.HHBlits,
                'binary': self.hhblits_binary_path,
                'format': 'a3m',
                'parser': parsers.parse_a3m},
        }

        unpairable_msas = []
        pairable_msas = []
        #for db_config in db_configs[0]:
            
        for db_config in self.db_configs:
            db_name = db_config['db_name']
            runner_type = db_config['runner']
            runner_config = runner_configs[runner_type]
            search_runner = runner_config['runner'](
                    binary_path=runner_config['binary'],
                    database_path=db_config['db_path'])

            msa_path = os.path.join(output_dir, db_name+'.'+runner_config['format'])

            raw_msa_result = pipeline.run_msa_tool(
                msa_runner=search_runner,
                input_fasta_path=fasta_path,
                msa_out_path=msa_path,
                msa_format=runner_config['format'],#a3m
                max_sto_sequences=db_config['max_hits'], #only for jackhmmer dbs
                use_precomputed_msas=True)

            raw_msa = raw_msa_result[runner_config['format']]
            parsed_msa = runner_config['parser'](raw_msa)

            if self.printout:
                print('MSA depth: {} - {} sequences'.format(db_name, len(parsed_msa)))

            if db_config['pairable'] == False:
                unpairable_msas.append(parsed_msa)
            else:
                pairable_msas.append(parsed_msa)        
        if self.easl_path != None:
            self.make_domtab(output_dir, fasta_path)
        #self.make_domtab(output_dir, fasta_path)

        return {'unpairable': unpairable_msas, 'pairable': pairable_msas}
    def get_MSAs_std(self):
        output_path = self.output_path
        #print(f'output path = {output_path}')
        for seq_name, seq_str in self.seq_dict.items():
            if self.printout:
                print(seq_name)
            output_dir = os.path.join(output_path, f'{seq_name}', 'std')
            if self.printout:
                print(output_dir)
            os.makedirs(output_dir, exist_ok=True)
            seq_msas = self.get_MSAs_for_seq_std(seq_str,
                                            output_dir)
            self.msa_dict[seq_name] = seq_msas
        return self.msa_dict

