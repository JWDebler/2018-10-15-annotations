import os
import re
import csv
from pathlib import Path

# This script expects a folder structure like this:
# 
# ├── candidates.py
# ├── output
# |     ├── 206-15
# |         ├── 206-15.deepsig.out
# |         ├── 206-15.effectorP.tsv
# |         ├── 206-15.genemark.proteins.fasta
# |         ├── 206-15.interproscan.tsv
# |         ├── 206-15.cazymes.txt <-- this one is optional
#
# scrape 'output' directory and use subfoldernames as sample names for further processing,
# this requires a folder structure of 'output/sampleID'
print("Reading content of output directory")
samples = [dI for dI in os.listdir('output') if os.path.isdir(os.path.join('output',dI))]
print("found following isolates:")
for sample in samples:
	print(sample)
print('============================================')

################ SETUP ################
#Conditions for effector candidates
effectorPscore = 0.8 # min EffectorP score
MWcutoff = 25  # molecular weight of mature (cleaved signal peptide) candidate
Cyscutoff = 2 # min number of cysteins
#######################################
for sample in samples:

	script_dir = Path.cwd()
	
	#Change this to where outputfiles are
	isolate_output_dir = script_dir / "output" / sample
	
	proteins = {}
	deepsig = {}
	signalP = {}
	signalPandDeepsig = {}

	effectorP = {}
	interproScan = {}
	cysteins = {}
	proteinMW = {}
	matureProtein = {}

	everything = {}
	candidates = {}
	cazymes = {}

	#Preprocessing
	print(sample,'- Processing all the raw data')

	#Read fasta file with protein sequences and IDs
	#And fill lists with default values
	input_file = sample+'.genemark.proteins.fasta'
	input_file_path = isolate_output_dir / input_file
	with open(input_file_path) as file:
		print(sample,'- reading proteins file')
		input = file.read().splitlines()
		for line in input:
			if not line.strip(): continue
			if line[0] == '>':
				name = line[1:]
				proteins[name] = ''
				signalP[name] = 0
				deepsig[name] = 0
				effectorP[name] = 0
				cazymes[name] = 'NA'

			else:
				proteins[name] += line

	print(sample,'- processed proteins data')

	#Read interproScan output
	input_file = sample+'.interproscan.tsv'
	input_file_path = isolate_output_dir / input_file
	with open(input_file_path) as file:
		print(sample,'- reading interproScan file')
		input = csv.reader(file, delimiter='\t')
		for line in input:
			id = line[0]
			if id in interproScan:
				interproScan[id].append(line[3:])
				if line[3] == 'SignalP_EUK':
					signalP[id] = int(line[7])
			else:
				interproScan[id] = [line[3:]]
				if line[3] == 'SignalP_EUK':
					signalP[id] = int(line[7])

	print(sample, '- processed interproScan data')

	#Read dbCAN2 CAZyme output
	input_file = sample+'.cazymes.txt'
	input_file_path = isolate_output_dir / input_file
	if input_file_path.is_file():
		with open(input_file_path) as file:
			print(sample,'- reading CAZyme file')
			input = csv.reader(file, delimiter='\t')
			next(input, None)
			for line in input:
				if int(line[5]) >= 2:
					cazymes[line[0]] = 'Yes'
				else:
					cazymes[line[0]] = 'No'
					
			print(sample,'- processed dbCAN cazyme data')
			
	else: 
		print('no cazymes file')
	
	

	#Read deepsig output
	input_file = sample+'.deepsig.out'
	input_file_path = isolate_output_dir / input_file
	with open(input_file_path) as file:
		print(sample,'- reading deepsig file')
		input = csv.reader(file, delimiter='\t')
		for line in input:
			if line[1] == 'SignalPeptide':
				deepsig[line[0]] = int(line[3])


	print(sample,'- processed deepsig data')

	#Read effectorP output
	input_file = sample+'.effectorP.tsv'
	input_file_path = isolate_output_dir / input_file
	with open(input_file_path) as file:
		print(sample,'- reading effectorP file')
		input = csv.reader(file, delimiter='\t')
		for line in input:
			if line[1] == 'Effector':
				effectorP[line[0]] = float(line[2])

	print(sample, '- processed effectorP data')

	#calculate cystein content and MW in kDA
	def cysteineContent(aaSequence):
		cystein = 0
		proteinMass = 0.0
		weights = dict(
			A = 71.03711,
			C = 103.00919,
			D = 115.02694,
			E = 129.04259,
			F = 147.06841,
			G = 57.02146,
			H = 137.05891,
			I = 113.08406,
			K = 128.09496,
			L = 113.08406,
			M = 131.04049,
			N = 114.04293,
			P = 97.05276,
			Q = 128.05858,
			R = 156.10111,
			S = 87.03203,
			T = 101.04768,
			V = 99.06841,
			W = 186.07931,
			Y = 163.06333)
		for aminoAcid in aaSequence:
			proteinMass += weights[aminoAcid]
			if aminoAcid == 'C':
				cystein += 1
		return cystein, round(proteinMass/1000,3)

	#Filling up 'everything' dictionary
	print(sample, '- Putting everything in one place')
	for e in proteins:
		everything[e] = ('deepsig', deepsig[e]), ('SignalP', signalP[e]), ('EffectorP', effectorP[e]),('Cazyme', cazymes[e]), ('Sequence', proteins[e])

	#Effector Prediction
	# Combine predictions where both SignalP and deepsig predict a signal peptide,
	# but if they can't agree on the cleavage site prefer SignalP over deepsig
	for key, value in signalP.items():
		if value > 0 and deepsig[key] > 0:
			if value == deepsig[key]:
				signalPandDeepsig[key] = value
			else:
				signalPandDeepsig[key] = value


	#Process cystein content and protein mass for candidates with signal peptides
	print(sample,'- calculating cystein content and MW for candidates')
	for key, value in signalPandDeepsig.items():
		mature = proteins[key][value:]
		temp = cysteineContent(mature)
		cysteins[key], proteinMW[key] = temp[0], temp[1]
		matureProtein[key] = mature

	print(sample,'- processed cystein content and MW for candidates')

	#defining ouput files
	output_name_candidates = sample+'.candidates.txt'
	output_file_candidates = isolate_output_dir / output_name_candidates

	output_name_candidates_mature_fasta = sample+'.candidates.mature.fasta'
	output_file_candidates_mature_fasta = isolate_output_dir / output_name_candidates_mature_fasta

	output_name_candidates_fasta = sample+'.candidates.fasta'
	output_file_candidates_fasta = isolate_output_dir / output_name_candidates_fasta

	#everything: SignalP[0/1], EffectorP[float], Cysteins[int], MW[float], aaSequence[string]
	for element in matureProtein:
		#If effectorP score > 0:
		if effectorP[element] >= effectorPscore:
			if proteinMW[element] <= MWcutoff:
				if cysteins[element] >= Cyscutoff:
					if element in interproScan:
						candidates[element] = effectorP[element], cysteins[element], proteinMW[element], cazymes[element], interproScan[element], matureProtein[element]
						print(element,'\t','SignalP=SP', '\t', 'EffectorP_Score:', effectorP[element],'\t', 'Cysteins:', cysteins[element], '\t', 'MW:', proteinMW[element], '\t', 'CAZyme:', cazymes[element], '\t', 'InterProScan_Results:', interproScan[element], file=open(output_file_candidates,"a"))
						print('>'+element, '\n'+matureProtein[element], file=open(output_file_candidates_mature_fasta,'a'))
						print('>'+element, '\n'+proteins[element], file=open(output_file_candidates_fasta,'a'))
					else:
						candidates[element] = effectorP[element], cysteins[element], proteinMW[element], cazymes[element], 'No InterProScan Hits', matureProtein[element]
						print(element,'\t','SignalP=SP', '\t', 'EffectorP_Score:', effectorP[element],'\t', 'Cysteins:', cysteins[element], '\t', 'MW:', proteinMW[element],  '\t', 'CAZyme:', cazymes[element], '\t', 'InterProScan_Results:', 'NONE', file=open(output_file_candidates,"a"))
						print('>'+element, '\n'+matureProtein[element], file=open(output_file_candidates_mature_fasta,'a'))
						print('>'+element, '\n'+proteins[element], file=open(output_file_candidates_fasta,'a'))

	print('============================================')
	print(sample,'============= Summary ============')
	print(sample,'- Selection of effector candidates:')
	print(sample,'- has Signal peptide (SignalP & deepsig)')
	print(sample,'- EffectorP score >=', effectorPscore)
	print(sample,'- Cysteins >=', Cyscutoff)
	print(sample,'- MW of mature protein <=', MWcutoff, 'kDa')
	print(sample,'- Effector candidates:',len(candidates))
	print('============================================')
	
