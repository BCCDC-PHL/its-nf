#!/usr/bin/env python3
import os, sys 
import pandas as pd
import argparse
from functools import partial
import re
from custom_html import HEAD, build_dbnote, build_table,build_table_title, build_row, FOOT, PAGEBREAK


EXPR_PRIORITY = re.compile(r'ATCC|LMG|type|NCTC')
ROW_NAMES = "subject_accession species bitscore percent_coverage percent_identity database_name extra_info".split()


def parse_db_csv(filepath):
	df = pd.read_csv(filepath, index_col=0)
	return df

def parse_fasta(fasta_path):
	fasta = []
	header = ''
	seq = ''
	with open(fasta_path, 'r') as f:
		for line in f:
			if line.startswith('>'):
				if header != '':
					fasta.append((header, seq))
				header = line.strip().lstrip('>')
				seq = ''
			else:
				seq += line.strip()
		fasta.append((header, seq))
	return fasta

def extract_descriptions(database_df):
	ncbi_path = os.path.join(database_df.loc['ncbi', 'PATH'], database_df.loc['ncbi','DBNAME'])

	headers, _ = zip(*parse_fasta(ncbi_path))

	def clean(string):
		string = string.replace("'",'')
		return string

	headers = [clean(x) for x in headers]

	def extract_info(string):
		string = " ".join(string.split()[3:])

		if re.search("strain|isolate", string):
			string = re.split("strain|isolate", string, maxsplit=1)[1]
		
		return re.split("(?:ITS|,)", string)[0]

	extra_info = [(x.split()[0], extract_info(x)) for x in headers ]

	df = pd.DataFrame(extra_info, columns=['subject_accession', 'extra_info']) 

	df = df.fillna('N/A')
	
	return df


def parse_blast(filepath):
	df = pd.read_csv(filepath)

	#df['database'] = os.path.basename(filepath).split('_')[1]

	df["species"] = df["species"].apply(lambda x: "" if x == "a" else x).astype(str)
	df = df.sort_values(['percent_coverage', 'percent_identity'], ascending=[False, False])

	df['percent_identity'] = df['percent_identity'].round(3)
	df['percent_coverage'] = df['percent_coverage'].round(3)


	df = df.drop_duplicates(subset=['query_seq_id', 'subject_accession', 'species', 'genus', 'percent_identity', 'percent_coverage', 'bitscore'])
	df = df[~df['species'].str.contains('uncultured')]

	return df.fillna('N/A')

def build_table_string(name, df, limit=20):
	build_row_partial = partial(build_row, row_names=ROW_NAMES)

	if df.shape[0] < limit:
		str_rows = df.apply(build_row_partial, axis=1)
		str_rows = '\n'.join(str_rows)
		str_table = build_table(name, ROW_NAMES, str_rows)

	else:
		N = df.shape[0]
		str_rows = df.iloc[0:limit].apply(build_row_partial, axis=1)
		str_rows = '\n'.join(str_rows)
		
		hidden_str_rows = df.iloc[limit+1:min(N, 300)].apply(build_row_partial, axis=1)
		hidden_str_rows = '\n'.join(hidden_str_rows)

		str_table = build_table(name, ROW_NAMES, str_rows, hidden_str_rows)
	
	return str_table

def main(args):
	local_blast_table = parse_blast(args.blast)
	ncbi_blast_table = parse_blast(args.ncbi)
	ncbi_blast_table['extra_info'] = ''

	database_df = parse_db_csv(args.db)
	DBNOTE = build_dbnote(database_df)

	extra_info = extract_descriptions(database_df)

	local_blast_table = local_blast_table.merge(extra_info, on='subject_accession', how='left')

	local_blast_dict = dict(list(local_blast_table.groupby('query_seq_id')))
	ncbi_blast_dict = dict(list(ncbi_blast_table.groupby('query_seq_id')))

	with open(args.output, "w") as outfile:
		outfile.write(HEAD)
		for name in set(local_blast_dict.keys()).union(ncbi_blast_dict.keys()):
			
			print(name)
			outfile.write(DBNOTE)

			outfile.write(build_table_title(name))

			if name in local_blast_dict:
				outfile.write(build_table_string(name, local_blast_dict[name]))

			if name in ncbi_blast_dict:
				outfile.write(build_table_string(name, ncbi_blast_dict[name]))
			
			outfile.write(PAGEBREAK)

		outfile.write(FOOT)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-b', '--blast', required=True, help='A single concatenated BLAST CSV table with hits from multiple samples and multiple database sources.')
	parser.add_argument('-n', '--ncbi', required=True, help='A single concatenated BLAST CSV table with hits from multiple samples from NCBI core_nt database.')
	parser.add_argument('-d', '--db', help='Database CSV file containing ID, DBNAME, and PATH columns.')
	parser.add_argument('-o', '--output', help='Output HTML report filename.')
	args = parser.parse_args()
	main(args)
