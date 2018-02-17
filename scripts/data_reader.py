"""
This class will read an excel file that contains somatic mutation data 
(sheet 1) and the genes (sheet 2). This will recreate the chart seen in
Kim et al.'s paper

Authors: Christian Choe, Min Cheol Kim
Create Date: 5/19/2016
"""

import collections
import os
import openpyxl

class data_reader:
	def __init__(self, file_name):
		#find the path to the data file
		self.file_name = file_name
		cur_path = os.path.dirname(__file__)
		new_path = os.path.relpath('../Data/' + file_name, cur_path)
		self.init_excel(new_path)

	#Opens the excel file and saves the two worksheets
	def init_excel(self, file_path):
		wb = openpyxl.load_workbook(file_path)
		excel_schart = wb.worksheets[0]
		excel_gchart = wb.worksheets[1]
		self.init_process_gchart(excel_gchart)
		self.init_process_schart(excel_schart)

	#Goes through the worksheets and gathers the important data
	def init_process_gchart(self, gchart):
		i = 2 #skip the first 2 rows because those are headers
		gene_list = {} #{key = coordinate; value = (gene_name, gene_frequency, prediction)}

		while gchart['A' + str(i + 1)].value:
			cur_row = gchart.rows[i]
			coordinate = cur_row[0].value.encode('ascii')
			mutation = '(' + coordinate[-3] + '->' + coordinate[-1] + ')'
			value = (cur_row[6].value.encode('ascii') + mutation,\
				cur_row[4].value.encode('ascii'), cur_row[5].value)
			gene_list[coordinate] = value
			i += 1
		self.gene_list = gene_list

	def init_process_schart(self, schart):
		i = 3 #skip the first 3 rows because those are headers
		mutation_list = collections.OrderedDict() #{key = gene_name; value = (mutation for each cell)}


		while schart['A' + str(i + 1)].value in self.gene_list.keys():
			cur_row = schart.rows[i]

			gene_name = self.gene_list.get(cur_row[0].value.encode('ascii'))[0]
			prediction = self.gene_list.get(cur_row[0].value.encode('ascii'))[2]
			temp_mutations = self.get_mutation_list(schart, cur_row, i + 1)

			if (gene_name not in mutation_list.keys()):
				#check if the prediction is not scored and if there is a next gene in the list
				if prediction == 'Not scored' and (schart['A' + str(i + 2)].value in self.gene_list.keys()):
					#check next cell below to see if it's the same gene since you don't want to store the 'Not scored' readings
					next_gene_name = self.gene_list.get(schart.rows[i + 1][0].value.encode('ascii'))[0]
					if next_gene_name != gene_name:
						mutation_list[gene_name] = temp_mutations #Gather list of somatic mutations for this gene
				else:
					mutation_list[gene_name] = temp_mutations #Gather list of somatic mutations for this gene
			i += 1
		self.mutation_list = mutation_list

	def get_mutation_list(self, schart, cur_row, row_num):
		temp_mutations = []
		wild_type = cur_row[1].value.encode('ascii')
		j = 3 #start of single-cell mutation data (skip the average cancer tissue)

		while schart.cell(row = row_num, column = j + 1).value:
			SNP = cur_row[j].value.encode('ascii')
			if SNP == wild_type:
				temp_mutations.append(0)
			elif SNP in 'ATCG':
				temp_mutations.append(2)
			elif SNP.isalpha():
				temp_mutations.append(1)
			else:
				temp_mutations.append('-') #equivalent to a '-' in the data file
			j += 1
		return temp_mutations

	def get_gene_mutations(self, gene_name):
		return self.mutation_list.get(gene_name)
