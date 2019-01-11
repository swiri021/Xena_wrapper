import xenaPython as xena
import pandas as pd
import numpy as np
import functools


## Decorator for Error Warning
class warning_decorator(object):
	def __init__(self,func_type=0):
		self.func_type = func_type

	def __call__(self, func):
		func_type =self.func_type

		@functools.wraps(func)
		def wrapper(*args, **kwargs):
			## ID mapper warning
			if func_type==0:
				try:
					#self.func(self, *args)
					values = func(*args, **kwargs)
					return values
				except:
					#raise
					raise ValueError("Cannot find EntrezIDs in ID set, Please insert correct EntrezIDs or other EntrezIDs")

			## Main function warning( Expr , CNV )
			elif func_type==1:
				values = func(*args, **kwargs)
				if values.empty:
					raise ValueError("There is no dataset with your list, Please use other EntrezIDs")
				else:
					return values

			## Main function warning( Mut )
			elif func_type==2:
				values = func(*args, **kwargs)
				if values==[]:
					raise ValueError("There is no dataset with your list, Please use other EntrezIDs")
				else:
					return values
		return wrapper

class load_TCGA:

	@warning_decorator(func_type=0)
	def id_mapper_to_entrez(self, df, input_list):
		#### ID mapping
		en_ids = pd.read_csv('TCGAlib/dataset/ID_table_homo_sapiens_20181211.csv').set_index('EnsemblID')[['EntrezID']]
		en_ids = en_ids[~en_ids.index.duplicated(keep='first')]
		crossed_ids = list(set(en_ids.index.tolist()).intersection(set(df.index.tolist())))
		en_ids = en_ids.loc[crossed_ids]

		id_map = pd.concat([en_ids, df], axis=1)
		id_map = id_map.reset_index().set_index('EntrezID')
		id_map = id_map[id_map.index.notnull()]
		id_map = id_map[~id_map.index.duplicated(keep='first')]
		id_map.index = id_map.index.astype(int).astype(str)

		return id_map.loc[input_list]

	@warning_decorator(func_type=1)
	def get_TCGA_expr(self, input_list):
		samples = xena.cohort_samples(self.host, [self.cohort], None) #### All Sample selection

		#### Create ID mapper
		xena_transcript = xena.dataset_field(self.host, self.expr_dataset) #### All Transcripts
		xena_transcript_df = pd.DataFrame(data=xena_transcript, columns =['EnsemblID'])
		xena_transcript_df['EnsemblID_edit'] = xena_transcript_df['EnsemblID'].apply(lambda x : x.split('.')[0])
		xena_transcript_df = xena_transcript_df.set_index('EnsemblID_edit')

		#### ID mapping & Get Ensembl ID from EntrezID
		id_map = self.id_mapper_to_entrez(xena_transcript_df, input_list)

		#### Retreive values from Xena
		values = xena.dataset_fetch(self.host, self.expr_dataset, samples, id_map['EnsemblID'].values.tolist()) # list of lists
		#xena_df = pd.DataFrame(data=values, index=input_list, columns=samples)

		xena_df = pd.DataFrame(data=values, index=input_list, columns=samples).T.astype(float)
		return xena_df

	@warning_decorator(func_type=1)
	def get_TCGA_cnv(self, input_list):

		samples = xena.cohort_samples(self.host, self.cohort, None) #### All Sample selection

		#### Create ID mapper
		pos_df = pd.read_csv('TCGAlib/dataset/position_info_hg38.csv', index_col=0)

		#### ID mapping
		id_map = self.id_mapper_to_entrez(pos_df, input_list)

		#### Get Ensembl ID from EntrezID
		id_map['chr'] = id_map['chr'].apply(lambda x:'chr'+x)
		id_map['start'].astype(int)
		id_map['end'].astype(int)

		#### Retreive values from Xena
		values = [xena.segmented_data_range(self.host, self.cnv_dataset, samples, id_map['chr'].loc[el], id_map['start'].loc[el], id_map['end'].loc[el]) for el in id_map.index.tolist()]
		cnv_arr = []
		for i,item in enumerate(input_list):
			temp_df = pd.DataFrame(data=values[i]['rows']['value'], columns=[item], index=values[i]['rows']['sampleID'])

			####There are duplicated samples in Xena, don't know the reason.
			temp_df = temp_df[~temp_df.index.duplicated(keep='first')].loc[samples]
			cnv_arr.append(temp_df)

		if len(input_list)==1:
			cnv_df = cnv_arr[0].astype(float)
		else:
			cnv_df = pd.concat(cnv_arr, axis=1).astype(float)

		return cnv_df

	@warning_decorator(func_type=2)
	def get_TCGA_mut(self, input_list):

		samples = xena.cohort_samples(self.host, self.cohort, None) #### All Sample selection

		#### Create ID mapper
		pos_df = pd.read_csv('TCGAlib/dataset/position_info_hg38.csv', index_col=0)

		#### ID mapping
		id_map = self.id_mapper_to_entrez(pos_df, input_list)

		#### Get Ensembl ID from EntrezID
		id_map['chr'] = id_map['chr'].apply(lambda x:'chr'+x)
		id_map['start'].astype(int)
		id_map['end'].astype(int)

		#### Retreive values from Xena
		values = [xena.sparse_data_range(self.host, self.mut_dataset, samples, id_map['chr'].loc[el], id_map['start'].loc[el], id_map['end'].loc[el]) for el in id_map.index.tolist()]
		result_arr = [pd.DataFrame(data=[values[i]['rows']['amino-acid'], values[i]['rows']['effect']], columns=values[i]['rows']['sampleID'], index=[item+'_AAC', item+'_Effect']).T for i,item in enumerate(input_list)]
		#mut_df = pd.concat(result_arr, axis=1)

		return result_arr

	def get_TCGA_surv(self):

		samples = xena.dataset_samples(self.host, self.surv_dataset, None)
		values = xena.dataset_fetch(self.host, self.surv_dataset, samples,['_EVENT','_TIME_TO_EVENT'])
		surv_df = pd.DataFrame(data=values, index=['Event', 'Time_to_event'], columns=samples).T
		return surv_df

	def __init__(self, host, user_cohort, prjID):
		self.host = host
		self.cohort=user_cohort

		self.cnv_dataset='%s/Xena_Matrices/%s.masked_cnv.tsv'%(prjID,prjID)
		self.surv_dataset = '%s/Xena_Matrices/%s.survival.tsv'%(prjID,prjID)
		self.expr_dataset = '%s/Xena_Matrices/%s.htseq_fpkm.tsv'%(prjID,prjID)
		self.mut_dataset = '%s/Xena_Matrices/%s.mutect2_snv.tsv'%(prjID,prjID)
		print self.host, self.cohort, prjID
