import fnmatch
import os
import argparse
import sys
import re
from ete3 import Tree
import time
from datetime import timedelta, date
import datetime
import csv
import json
from statistics import mode




def nwkToAnnotNameNWK(t, j_d, trait, inMeta_name, pangolin, indexToGene, geneToIndex):
	'''
	Input: ete3 tree with node names that have 'trait' of clade specified in j_d
	Outputs: tree with trait appened to node names
	'''
	firstDate = date.fromisoformat('2030-01-01')
	lastDate =  date.fromisoformat('1900-01-01')

	sampPangolin_d = {}
	sampDate_d = {}
	strain_index = "NA"
	pangolin_index = "NA"
	inMeta_open = open(inMeta_name, "r")
	for line in inMeta_open:
		line_l = re.split('\t|,', line.strip())
		if strain_index != "NA":
			tipName = line_l[strain_index]
			tipDate = line_l[date_index]

			sampDate_d[tipName] = tipDate 

			if tipDate != "":
				isoDate = date.fromisoformat(tipDate)

				if isoDate < firstDate:
					firstDate = isoDate
				if isoDate > lastDate:
					lastDate = isoDate

			if pangolin == 'metadata':
				sampPangolin = line_l[pangolin_index]
				sampPangolin_d[tipName] = sampPangolin

		elif "strain" in line_l and "date" in line_l:
			strain_index = line_l.index("strain")
			date_index = line_l.index("date")

		elif "Name" in line_l and "Date" in line_l:
			strain_index = line_l.index("Name")
			date_index = line_l.index("Date")

		else:
			sys.exit(line_l)
			#"Must have 'strain' and 'date' in first line of tab seperated"

		if pangolin == 'metadata' and pangolin_index == "NA":
			if "pangolin_lin" in line_l:
				pangolin_index = line_l.index("pangolin_lin")
			if "lineage" in line_l:
				pangolin_index = line_l.index("lineage")
	
	if pangolin != 'metadata': # pangolin output lineage report file name
		inPangolin_open = open(pangolin, "r")
		for line in inPangolin_open:
			line_l = line.split(",")
			if line_l[2] != "probability": #is not first line
				if float(line_l[2]) > 0.9:
					sampPangolin = line_l[1]
				else:
					print("Prob does not meet threshold", tipName, float(line_l[2]) )
					sampPangolin = 'NA'
				tipName = line_l[0]
				sampPangolin_d[tipName] = sampPangolin


	for node in t.traverse("preorder"):
		if "NODE_" not in node.name:
			if node.name in sampPangolin_d:
				pID = sampPangolin_d[node.name]
			else: 
				pID = "NA"
			temp = "_" + pID + "_" + sampDate_d[node.name] + "_"
		else:
			temp = "_"

		if trait == "aa_muts":
			#for gene in geneToIndex:
			for index in range(12):
				gene = indexToGene[index]
				if gene in j_d['nodes'][node.name][trait]:
					for mut in j_d['nodes'][node.name][trait][gene]:
						usable = mut.replace(".", "J").replace("-", "X").replace("*", "J") #gaps represented by X
						temp += usable+"."
				temp += "-"
		else:
			temp += str(j_d['nodes'][node.name][trait])
			#node_annot = temp
		node.name = node.name + temp
	return(t, firstDate, lastDate)


def assignToSpecAA(t, mutList, geneToIndex, logNotes_open):
	'''
	(assignment_d: key: leaf node name; value: clade) and heierarchy (heiarchy_d: key:child clade; value:parent clade )
	'''
	clade_s = set()
	assignment_d = {}
	heiarchy_d = {}
	cladesAssinged_d ={}

	currentClade = 'anc'
	heiarchy_d['anc'] = 'NA'
	clade_s.add('anc')

	for node in t.traverse("preorder"):
		assignment_d[node.name] = currentClade

	for node in t.traverse("preorder"):

		node_name = node.name

		nodeMut_l = node_name.split("_")[-1].split("-")
		currentClade = assignment_d[node_name]

		mutOfInterst_count = 0

		for gene_and_mutToFind in mutList:
			for gene in geneToIndex:
				if gene == gene_and_mutToFind[:len(gene)]:
					mutToFind = gene_and_mutToFind[len(gene):]

					nodeMut_gene_catString = nodeMut_l[geneToIndex[gene]]
					nodeMut_gene_l = nodeMut_gene_catString.split(".")
					if nodeMut_gene_l != ['']:
						for nodeMut_one in nodeMut_gene_l:
							if nodeMut_one != '':

								if mutToFind[0] == "*" or mutToFind[0] == nodeMut_one[0]:
									if mutToFind[-1] == "*" or mutToFind[-1] == nodeMut_one[-1]:
										if mutToFind[1:-1] == nodeMut_one[1:-1]:
											mutOfInterst_count += 1
											mutFound = nodeMut_one
		if mutOfInterst_count == 1:
			if mutFound not in cladesAssinged_d:
				cladesAssinged_d[mutFound] = 0
			cladesAssinged_d[mutFound] += 1

			newClade = mutFound+ "_" + str(cladesAssinged_d[mutFound])
			clade_s.add(newClade)

			heiarchy_d[newClade] = currentClade
			assignment_d[node_name] = newClade

			for des_node in node.iter_descendants("postorder"):
				assignment_d[des_node.name] = newClade

		if mutOfInterst_count > 1:
			#print("found more than 1 mut of interst", node_name)
			logNotes_open.write("\n"+"found more than 1 mut of interst"+"\n")
			logNotes_open.write(node_name+"\n")



	return(assignment_d, heiarchy_d, clade_s)


def assignToTraits(t, ofInterst_l = []):
	'''
	(assignment_d: key: leaf node name; value: clade) and heierarchy (heiarchy_d: key:child clade; value:parent clade )
	calde naming: trait_index
	'''
	clade_s = set()
	assignment_d = {}
	heiarchy_d = {}
	cladesAssinged_d ={}

	currentClade = 'anc'
	heiarchy_d['anc'] = 'NA'
	clade_s.add('anc')

	for node in t.traverse("preorder"):
		assignment_d[node.name] = currentClade

	for node in t.traverse("preorder"):

		if not node.is_root():

			node_name = node.name

			traitOfNode = node_name.split("_")[-1]
			currentClade = assignment_d[node_name]


			#define new calde because does not match parent
			if currentClade.split("_")[0] != traitOfNode:
				if traitOfNode not in cladesAssinged_d:
					cladesAssinged_d[traitOfNode] = 0
				cladesAssinged_d[traitOfNode] += 1
				newClade = traitOfNode + "_" + str(cladesAssinged_d[traitOfNode])
				clade_s.add(newClade)

				heiarchy_d[newClade] = currentClade
				assignment_d[node_name] = newClade

				for des_node in node.iter_descendants("postorder"):
					assignment_d[des_node.name] = newClade


	return(assignment_d, heiarchy_d, clade_s)




def assignCladeToLin(assignment_old_d, heiarchy_old_d, clade_old_s):


	reassignCladeNames ={}
	clade_s = set()
	assignment_d = {}
	heiarchy_d = {}
	
	for clade in clade_old_s:
		tipLinForClade = []
		for node in assignment_old_d:
			if 'NODE' not in node:
				if assignment_old_d[node] == clade:
					#print(node)
					tipLinForClade = tipLinForClade + [node.split("_")[-3]]
		if clade != 'anc':
			try:
				lineage = mode(tipLinForClade)
			except ValueError:
				lineage = "2lin"
			newClade = lineage + "_" + clade
		else:
			newClade = clade

		
		reassignCladeNames[clade] = newClade
		clade_s.add(newClade)

	for child in heiarchy_old_d:
		parent = heiarchy_old_d[child]
		if parent != 'NA':
			heiarchy_d[reassignCladeNames[child]] = reassignCladeNames[parent]
		else:
			heiarchy_d[reassignCladeNames[child]] = parent
	for samp in assignment_old_d:
		assignment_d[samp] = reassignCladeNames[assignment_old_d[samp]]

	return (assignment_d, heiarchy_d, clade_s)


def countAbudanceFromNames_byWeek(assignment_d, clade_s, startDate, endDate, delta, tipLog_name):
	abundances_d = {} # key: week; value: dict of key:clade; value: count
	weekToDate_d = {}
	#assignment_d: key: node name; value: clade

	week_l = []

	tipLog_open = open(tipLog_name, "w")
	outLine = "	".join(["Week", "Clade", "sample_withAnnot", "sample_inputID"]) + "\n"
	tipLog_open.write(outLine)

	currentStart = startDate
	currentEnd = currentStart + delta
	week = 0
	weekName = str(week)
	week_l = week_l + [weekName]
	abundances_d[weekName] = {}

	while endDate >= currentEnd:

		lastWeekName = weekName
		week += 1
		weekName = str(week)
		week_l = week_l + [weekName]

		weekToDate_d[weekName] = currentStart

		abundances_d[weekName] = {} 

		for clade in abundances_d[lastWeekName]:
			abundances_d[weekName][clade] = 1


		for tip in assignment_d.keys():
			if "NODE_" not in tip and "Wuhan" not in tip:
				tip_date = date.fromisoformat(tip.split("_")[-2])

				if tip_date < currentEnd and tip_date >= currentStart:

					clade = assignment_d[tip]
					if clade not in abundances_d[weekName]:
						abundances_d[weekName][clade] = 1
					abundances_d[weekName][clade] += 1
					if clade != 'anc':
						tipNoAnnot = ""
						for name in tip.split("_")[:-3]: 
							tipNoAnnot = tipNoAnnot + "_"  + name
						outLine = "	".join([ weekName, clade, tip, tipNoAnnot[1:]]) + "\n"
						tipLog_open.write(outLine)


		currentEnd += delta
		currentStart += delta
	tipLog_open.close()

	#print(list(clade_s))
	noFurtherAbudance = list(clade_s.copy())
	for weekName in (reversed(week_l)):

		noFurtherAbudance_last = noFurtherAbudance.copy()
		for clade in noFurtherAbudance_last:
			if clade in abundances_d[weekName]:

				if abundances_d[weekName][clade] == 1 and clade != 'anc':
					abundances_d[weekName][clade] = 0
				else:
					noFurtherAbudance.remove(clade)

	return(abundances_d, weekToDate_d)


def main():



	##########################  parse user arguments
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-n', '--inNextstrain', required=True, type=str, help="nextstrain results with tree.nwk and [traitOfInterst].json")
	parser.add_argument('-m', '--inMeta', required=True, type=str, help="metadata tsv with 'strain'	and 'date'cols, optional: cols of trait of interst; and pangolin col named: 'lineage' or 'pangolin_lin'")
	parser.add_argument('-p', '--inPangolin', required=False, type=str, default = "metadata", help="pangolin output lineage_report.csv file, if argument not supplied looks in inMeta for col with 'pangolin_lin' or 'lineage'")

	parser.add_argument('-f', '--traitOfInterstFile', required=False, type=str, default="aa_muts.json",  help="name of nextstrain [traitOfInterst].json in 'inNextstrain' folder")
	parser.add_argument('-k', '--traitOfInterstKey', required=False, type=str, default="aa_muts",  help="key for trait of interst in json file")
	parser.add_argument('-aa', '--aaVOClist', required=False, nargs='+', help="list of aa of interest in form [GENE][*ORAncAA][site][*ORtoAA] ex. S*501*, gaps represed by X")



	parser.add_argument('-oDir', '--outDirectory', required=False, default ="./", type=str, help="folder for output")
	parser.add_argument('-oP', '--outPrefix', required=True, type=str, help="prefix of out files withen outDirectory")


	parser.add_argument('-t', '--timeWindow', required=False, type=str, default="7", help="number of days for sampling window")
	parser.add_argument('-s', '--startDate', required=False, type=str, default='2020-03-01', help="start date in iso format YYYY-MM-DD or 'firstDate' which is in metadata")
	parser.add_argument('-e', '--endDate', required=False, type=str, default='lastDate', help="end date in iso format YYYY-MM-DD or 'lastDate' which is in metadata")

	#parser.add_argument('-w', '--MINTIME', required=False, type=str, default="30", help="minimum time point to start plotting")
	#parser.add_argument('-min', '--MINTOTALCOUNT', required=False, type=str, default="10", help="minimum total count for group to be included")

	

	args = parser.parse_args()

	##########################  file names

	inTree_name = os.path.join(args.inNextstrain, "tree.nwk")
	
	
	inJSON_name = os.path.join(args.inNextstrain, args.traitOfInterstFile)
	inMeta_name = args.inMeta
	# if args.inPangolin is None:
	# 	pangolin = 'metadata'
	# else:
	pangolin = args.inPangolin
	
	outCladeHierarchy_name = os.path.join(args.outDirectory, args.outPrefix +"_defineCountClades", "cladeHierarchy.csv")
	outCounts_name = os.path.join(args.outDirectory, args.outPrefix +"_defineCountClades", "abundances.csv")
	tipLog_name = os.path.join(args.outDirectory, args.outPrefix +"_defineCountClades", "tipsWithMut.txt")

	logNotes_name = os.path.join(args.outDirectory, args.outPrefix + "_defineCountClades", "log.txt")


	# weekTreeFile_folder = os.path.join(args.outDirectory, args.outPrefix+ "_PrunedTrees")
	# weekTreeFile_prefix = os.path.join(weekTreeFile_folder, "PrunedWeek_")



	if not os.path.exists(args.outDirectory):
		os.makedirs(args.outDirectory)
	# if not os.path.exists(weekTreeFile_folder):
	#  	os.makedirs(weekTreeFile_folder)
	dataOutDir = os.path.join(args.outDirectory, args.outPrefix + "_defineCountClades")
	if not os.path.exists(dataOutDir):
	 	os.makedirs(dataOutDir)


	######################### params to referance

	geneToIndex = {"E":0, "M":1, "N":2, "ORF1a":3, "ORF1b":4, "ORF3a":5, "ORF6":6, "ORF7a":7, "ORF7b":8, "ORF8":9, "ORF9b":10, "S":11}
	indexToGene = {0:"E", 1:"M", 2:"N", 3:"ORF1a", 4:"ORF1b", 5:"ORF3a", 6:"ORF6", 7:"ORF7a", 8:"ORF7b", 9:"ORF8", 10:"ORF9b", 11:"S"}


	delta= timedelta(days=int(args.timeWindow))


	#########################  load nwk files

	print(inTree_name)
	if not os.path.exists(inTree_name):
		sys.exit("missing input tree")
	t = Tree(inTree_name, format = 3)


	j_d = json.load(open(inJSON_name))

	
	t, firstDate, lastDate = nwkToAnnotNameNWK(t, j_d, args.traitOfInterstKey, inMeta_name, pangolin, indexToGene, geneToIndex)

	if args.startDate == "firstDate":
		startDate = firstDate
	else:
		startDate =  date.fromisoformat(args.startDate)
	if args.endDate == "lastDate":
		endDate = lastDate
	else:
		endDate =  date.fromisoformat(args.endDate)	



	######################### make clades assignments (assignment_d: key: leaf node name; value: clade) and heierarchy (heiarchy_d: key:child clade; value:parent clade )

	logNotes_open = open(logNotes_name, "w")
	if args.traitOfInterstKey == 'aa_muts':
		if args.aaVOClist is None:
			#aa_mut = ['S*501*', 'ORF1b*1183*', 'S*484*', 'S*452*', 'S*681*', 'S*13', 'S*152']
			#aa_mut = ['S*484*', 'S*501*'] #default sites of interst
			aa_mut = ["S*484*", "S*501*", "S*13*", "S*452*", "S*477*", "N*205*", "S*253*"]
			# aa_mut = ["S*484*", "S*501*", "S*13*", "S*452*", "S*477*", "S*253*"]
			# aa_mut = ['S*484*']

		else:
			aa_mut = args.aaVOClist

		print("\n")
		print("Searching for muations with:", aa_mut)

		
		logNotes_open.write("\n"+"Searching for muations with:"+str(aa_mut)+"\n")


		assignment_d, heiarchy_d, clade_s = assignToSpecAA(t, aa_mut, geneToIndex, logNotes_open)
	else:
		logNotes_open.write("\n"+"Finding clades defined by"+str(args.traitOfInterstKey)+"\n")
		assignment_d, heiarchy_d, clade_s = assignToTraits(t)

	assignment_d, heiarchy_d, clade_s = assignCladeToLin(assignment_d, heiarchy_d, clade_s)


	######################### count abundances based on assignment_d


	abundances_d, weekToDate_d = countAbudanceFromNames_byWeek(assignment_d, clade_s, startDate, endDate, delta, tipLog_name)

	#abundances_d, heiarchy_d = clipPreParent(abundances_d, heiarchy_d)


	# ##########################  write abundances over time

	outCounts_open = open(outCounts_name, "w")
	headLine = "names,times,abundances,date\n"
	outCounts_open.write(headLine)

	for week in abundances_d:
		for name in abundances_d[week]:
			outLine = ",".join([name, str(week.replace("week", "")), str(abundances_d[week][name]), str(weekToDate_d[week])]) + "\n"
			outCounts_open.write(outLine)
	outCounts_open.close()

	##########################  write CladeHierarchy to file

	outCladeHierarchy_open = open(outCladeHierarchy_name, "w")
	headLine = "names,parents \n"
	outCladeHierarchy_open.write(headLine)
	for child in heiarchy_d:
		outLine = child + "," + heiarchy_d[child] + "\n"
		outCladeHierarchy_open.write(outLine)

	outCladeHierarchy_open.close()


	#print("FINAL heiarchy_d")
	#print(heiarchy_d)
	logNotes_open.write("\n"+"Final heiarchy_d"+"\n")
	logNotes_open.write(str(heiarchy_d)+"\n")


if __name__ == "__main__":
	MIN_PYTHON = (3, 7)
	if sys.version_info < MIN_PYTHON:
		sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)

	print("Loaded imports")
	main()