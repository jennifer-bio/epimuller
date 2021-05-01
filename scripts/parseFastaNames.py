import fnmatch
import os
import argparse
import sys



def main():

	# parse user arguments
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-iF', '--inFasta', required=True, type=str, help="full metadata file with ISO date after last '|' in name")
	parser.add_argument('-oM', '--outMeta', required=True, type=str, help="output for metadata file")
	parser.add_argument('-oF', '--outFasta', required=True, type=str, help="output for fasta file")
	parser.add_argument('-p', '--inPangolin', required=False, type=str, default = "metadata", help="pangolin output lineage_report.csv file, if argument not supplied adds ? in 'lineage' col")


	args = parser.parse_args()

	inFasta_open = open(args.inFasta, "r")
	outFasta_open = open(args.outFasta, "w")
	outMeta_open = open(args.outMeta, "w")

	lineage_d = {}
	if args.inPangolin is not None:
		haveLineage = True
		inPangolin_open = open(args.inPangolin, "r")
		for line in inPangolin_open:
			line_l = line.strip().split(",")
			name = line_l[0].replace("|", "_")
			lineage = line_l[1]
			lineage_d[name] = lineage
	else:
		haveLineage = False


	outMetaLine = "strain	virus	date	region	country_exposure	date_submitted	lineage\n"
	outMeta_open.write(outMetaLine)

	inName = 'NA'

	for line in inFasta_open:
		if ">" in line:
			if inName != 'NA' and 'Wuhan' not in inName:
				outFasta_open.write(outFastaLine)
				outFasta_open.write(seq)
				outMeta_open.write(outMetaLine)

			inName = line.strip()[1:]

			date = inName.split("|")[-1]
			outName = inName.replace("|", "_")
			outFastaLine = ">" + outName + "\n"
			if haveLineage:
				outMetaLine = "	".join([outName, "ncov", date, "North America", "?	2030-03-22", lineage_d[outName] +"\n"]) 
			else:
				outMetaLine = "	".join([outName, "ncov", date, "North America", "?	2030-03-22	?\n"]) 
			
			seq = ""
		else:
			seq += line
	outFasta_open.write(outFastaLine)
	outFasta_open.write(seq)
	outMeta_open.write(outMetaLine)



if __name__ == "__main__":
	main()


