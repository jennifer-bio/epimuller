import fnmatch
import os
import argparse
import sys
import re
from ete3 import Tree
import time
from datetime import timedelta, date
import csv
import json
from statistics import mode
import svgwrite #conda install svgwrite
import random
import cairo
import cairosvg



def main():


	MIN_PYTHON = (3, 7)
	if sys.version_info < MIN_PYTHON:
		sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)

	print("Loaded imports")

	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)


	##########################  parse user arguments
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	global_group  = parser.add_argument_group('Options for full repot')

	global_group.add_argument('-oDir', '--outDirectory', required=False, default ="./", type=str, help="folder for output")
	global_group.add_argument('-oP', '--outPrefix', required=True, type=str, help="prefix of out files withen outDirectory")


	defClades_group  = parser.add_argument_group('Options passed to epimuller-define')


	treeAndtraits = parser.add_mutually_exclusive_group(required=True)
	treeAndtraits.add_argument('-n', '--inNextstrain', type=str, help="nextstrain results with tree.nwk and [traitOfInterst].json")
	treeAndtraits.add_argument('-a', '--annotatedTree', type=str, help="nexus file name  and [traitOfInterst].json")

	#defClades_group.add_argument('-n', '--inNextstrain', required=True, type=str, help="nextstrain results with tree.nwk and [traitOfInterst].json")
	defClades_group.add_argument('-m', '--inMeta', required=True, type=str, help="metadata tsv with 'strain'	and 'date'cols, optional: cols of trait of interst; and pangolin col named: 'lineage' or 'pangolin_lin'")
	defClades_group.add_argument('-p', '--inPangolin', required=False, type=str, default = "metadata", help="pangolin output lineage_report.csv file, if argument not supplied looks in inMeta for col with 'pangolin_lin' or 'lineage'")
	defClades_group.add_argument("--noPangolin", action="store_true", help="do not add lineage to cade names")

	defClades_group.add_argument('-f', '--traitOfInterstFile', required=False, type=str, default="aa_muts.json",  help="name of nextstrain [traitOfInterst].json in 'inNextstrain' folder")
	defClades_group.add_argument('-g', '--geneBoundry', required=False, type=str, help="json formated file specifing start end postions of genes in alnment for annotatedTree with aa_muts option")

	defClades_group.add_argument('-k', '--traitOfInterstKey', required=False, type=str, default="aa_muts",  help="key for trait of interst in json file or annotated tree file for aa with 'mutations' annotation, use 'aa_muts'")
	defClades_group.add_argument('-mut', '--VOClist', required=False, nargs='+', help="list of aa of interest in form [GENE][*ORAncAA][site][*ORtoAA] ex. S*501*, gaps represed by X")


	defClades_group.add_argument('-t', '--timeWindow', required=False, type=str, default="7", help="number of days for sampling window")
	defClades_group.add_argument('-s', '--startDate', required=False, type=str, default='2020-03-01', help="start date in iso format YYYY-MM-DD or 'firstDate' which sets start date to first date in metadata")
	defClades_group.add_argument('-e', '--endDate', required=False, type=str, default='lastDate', help="end date in iso format YYYY-MM-DD or 'lastDate' which sets end date as last date in metadata")

	

	##########################  command line args for drawMuller.py

	drawing_group  = parser.add_argument_group('Options passed to epimuller-draw')

	drawing_group.add_argument('-mt', '--MINTIME', required=False, type=str, default="30", help="minimum time point to start plotting")
	drawing_group.add_argument('-min', '--MINTOTALCOUNT', required=False, type=str, default="50", help="minimum total count for group to be included")

	drawing_group.add_argument('-c', '--cases_name', required=False, type=str, help="file with cases - formated with 'date' in ISO format and 'confirmed_rolling' cases, in tsv format")

	drawing_group.add_argument('-l', '--xlabel', required=False, type=str, choices = ["date", "time"], default="date", help="Format of x axis label: ISO date format or timepoints from start")
	drawing_group.add_argument('-lp', '--labelPosition', required=False, type=str, default="Right", choices = ["Right", "Max", "Start", "End"], help="choose position of clade labels")

	drawing_group_page  = parser.add_argument_group('Options passed to epimuller-draw for page setup')
	drawing_group_page.add_argument('--WIDTH', required=False, type=str, default="1500", help="WIDTH of page (px)")
	drawing_group_page.add_argument('--HEIGHT', required=False, type=str, default="1000", help="HEIGHT of page (px)")
	drawing_group_page.add_argument('--LEGENDWIDTH', required=False, type=str, default="220", help="LEGENDWIDTH to the right of plotting area (px)")
	drawing_group_page.add_argument('--MARGIN', required=False, type=str, default="60", help="MARGIN around all sides of plotting area (px)")
	drawing_group_page.add_argument('--FONTSIZE', required=False, type=str, default="26")
	drawing_group_page.add_argument('--LABELSHIFT', required=False, type=str, default="15", help="nudge label over by LABELSHIFT (px)")
	


	args = parser.parse_args()


	#call with script
	# commandCallDefine = "python ../../epiMuller/scripts/defineAndCountClades.py"
	# commandCallDraw = "python ../../epiMuller/scripts/drawMuller.py"


	#call with entry_points

	commandCallDefine = "epimuller-define"
	commandCallDraw = "epimuller-draw"


	########################### call defineAndCountClades.py

	if args.inNextstrain is not None:
		treeAndtraits = "--inNextstrain " + args.inNextstrain
	else:
		treeAndtraits = "--annotatedTree " + args.annotatedTree

	if args.noPangolin:
		noPangolin = " --noPangolin "
	else:
		noPangolin = ""

	if args.geneBoundry is not None:
		geneBoundry = " --geneBoundry " + args.geneBoundry + " "
	else:
		geneBoundry = ""


	oscommand = " ".join([commandCallDefine, "--outDirectory", args.outDirectory , "--outPrefix", args.outPrefix ,
	 treeAndtraits , "--inMeta", args.inMeta , "--inPangolin", args.inPangolin , 
	 "--traitOfInterstFile", args.traitOfInterstFile , "--traitOfInterstKey", args.traitOfInterstKey , geneBoundry + noPangolin + "--timeWindow", args.timeWindow ,
	 "--startDate", args.startDate ,"--endDate", args.endDate])


	if args.VOClist is not None:
		oscommand = oscommand + " --VOClist"
		for aa in args.VOClist:
			oscommand += " " + aa


	print("\n ### Call command ###")
	print(oscommand)
	print("\n")
	if os.system(oscommand) != 0:
		sys.exit("epimuller-define failed")


	##########################  call drawMuller.py


	outCladeHierarchy_name = os.path.join(args.outDirectory, args.outPrefix +"_defineCountClades", "cladeHierarchy.csv")
	outCounts_name = os.path.join(args.outDirectory, args.outPrefix +"_defineCountClades", "abundances.csv")
	plot_folder = os.path.join(args.outDirectory, args.outPrefix+ "_PLOTS", "")


	# oscommand = " ".join(["python 00_scripts/drawMuller.py --parentHierarchy_name", outCladeHierarchy_name, "--abundance_name", outCounts_name, "--outFolder", plot_folder,
	# 	"--xlabel", args.xlabel, "--labelPosition", args.labelPosition, "--MINTIME", args.MINTIME, "--MINTOTALCOUNT", 
	# 	args.MINTOTALCOUNT, "--xlabel", args.xlabel , "--labelPosition", args.labelPosition])
	oscommand = " ".join([commandCallDraw, "--parentHierarchy_name", outCladeHierarchy_name, "--abundance_name", outCounts_name, "--outFolder", plot_folder,
		"--xlabel", args.xlabel, "--labelPosition", args.labelPosition, "--MINTIME", args.MINTIME, "--MINTOTALCOUNT", 
		args.MINTOTALCOUNT, "--xlabel", args.xlabel , "--labelPosition", args.labelPosition, 
		"--WIDTH", args.WIDTH, "--HEIGHT", args.HEIGHT, "--LEGENDWIDTH", args.LEGENDWIDTH, "--MARGIN", args.MARGIN , "--FONTSIZE", args.FONTSIZE, "--LABELSHIFT", args.LABELSHIFT])

	if args.cases_name is not None:
		oscommand = oscommand + " --cases_name " + args.cases_name
	
	print("\n ### Call command ###")
	print(oscommand)
	print("\n")
	if os.system(oscommand) != 0:
		sys.exit("epimuller-draw failed")
	


if __name__ == "__main__":
	MIN_PYTHON = (3, 7)
	if sys.version_info < MIN_PYTHON:
		sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)

	print("Loaded imports")
	main()