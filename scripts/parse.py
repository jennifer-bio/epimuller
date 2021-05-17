
import fnmatch
import os
import argparse
import sys
import re

from treeClass import Node, parseTree

def getTreeRawString(tree_fileName, treeformat):
	'''Extracts string that defines free from newick or nexus file'''
	treeformat_d = parseFormat(treeformat)


	tree_fileOpen = open(tree_fileName, "r")

	if treeformat_d["fileFormat"] == "nexus":
		nextLine = False

		for line in tree_fileOpen:
			if 'Begin Trees;' in line:
				nextLine = True
			elif nextLine:
				treeRawString = "=".join(line.strip().split("=")[1:])
				nextLine = False
	elif treeformat_d["fileFormat"] == "newick":
		treeRawString = ""
		for line in tree_fileOpen:
			treeRawString += line.strip()
	else:
		"do not understand tree fileFormat"
		sys.exit()
 
	tree_fileOpen.close()

	if treeRawString[-1] == ";":
		treeRawString = treeRawString[:-1]

	return treeRawString

def matchParenCommas(treeString):
	'''
	input string
	output dictionaries of openParen_d: key:index of open parentheses, value: index of close parentheses; 
						   closeParn_d  key:index of close parentheses, value: index of open parentheses; 
						   commaParn_d:  key:index of open parentheses, value: list of indices of commas in that parenthesis levels
	'''
	stack_paren = []
	stack_comma = []
	openParen_d = {}
	closeParn_d = {}
	commaParn_d = {}
	notInAnnot = True

	for i, c in enumerate(treeString):
		if c == "(":
			stack_paren.append(i)
			stack_comma.append([])
		elif c == "," and notInAnnot:
			stack_comma[-1].append(i)
		elif c == ")":
			openParen = stack_paren.pop()
			openParen_d[openParen] = i
			closeParn_d[i] = openParen

			commaList = stack_comma.pop()
			commaParn_d[openParen] = commaList
		elif c == "[" or  c == "]" or c == "{" or c == "}":
			notInAnnot = not notInAnnot

	return(openParen_d, closeParn_d, commaParn_d)



def parseFormat(treeformat):
	'''
	checks for proporties of tree format list, that deviate from default; and returns dictionary
	'''
	
	treeformat_d = {"fileFormat":"newick", "branchLen":True, "internalNames":True, "leafNames":True, "annotation":False, "annotationType":None}

	if "nexus" in treeformat:
		treeformat_d["fileFormat"] = "nexus"


	if "noBranchLen" in treeformat:
		treeformat_d["branchLen"] = False


	if "treetimeAnnot" in treeformat:
		treeformat_d["annotation"] = True
		treeformat_d["annotationType"] = "treetimeAnnot"
	elif "hyphyAnnot" in treeformat:
		treeformat_d["annotation"] = True
		treeformat_d["annotationType"] = "hyphyAnnot"

	return (treeformat_d)


def buildTree(treeString, treeformat, openParen_d, closeParn_d, commaParn_d, subTstart, subTend, parent):
	'''
	builds tree by nodes pointing to parent and child based on input tree tree string
	returns root Node object
	'''


	currentNode = Node(parent, None, None)
	currntTreeString = treeString[subTstart+1:subTend]


	if "(" in currntTreeString: 
		#internal nodes

		firstParen = min(openParen_d.keys())
		closeParen = openParen_d.pop(firstParen)
		slitCommas = commaParn_d[firstParen]

		currentNodeString = treeString[closeParen+1:subTend]
		
		allSplits = [firstParen] + slitCommas + [closeParen]
		for i in range(len(allSplits)-1):
			childNode = buildTree(treeString, treeformat, openParen_d, closeParn_d, commaParn_d, allSplits[i], allSplits[i+1], currentNode)
			currentNode.children.append(childNode)

	else:
		#leaf
		currentNodeString = currntTreeString


	# parse currentNodeString name for details
	treeformat_d = parseFormat(treeformat)

	if treeformat_d["branchLen"]:
		if ":" in currentNodeString:
			currentNodeLen = currentNodeString.split(":")[1]		
			currentNodeName = currentNodeString.split(":")[0]
		else: #should only be in root for unrooted tree
			currentNodeLen = 0
			currentNodeName = currentNodeString
	else:
		currentNodeLen = 1
		currentNodeName = currentNodeString

	if treeformat_d["annotation"]:
		if treeformat_d["annotationType"] == "hyphyAnnot":
			annotation_l = currentNodeName.split("{")[-1].split(",")
			setattr(currentNode, "labels_l", annotation_l)
			currentNode.annotKeys.append("labels_l")
			currentNodeName = currentNodeName.split("{")[0]

		elif treeformat_d["annotationType"] == "treetimeAnnot":
			if "[" in currentNodeLen: #note: no annotation assigned to root 
				annotation_str = currentNodeLen.split("[")[1].replace("&", "")
				annotKey = ""
				annotValue = ""
				addToKey = True
				outOfQuote = True
				annotation_d = {}
				for c in annotation_str:
					if c == '"':
						outOfQuote = not outOfQuote

					if c == "=":
						addToKey = False
					elif (c == "," and outOfQuote) or c == "]":
						addToKey = True
						setattr(currentNode, annotKey, annotValue.replace('"', ''))
						currentNode.annotKeys.append(annotKey)
						annotKey = ""
						annotValue = ""
					elif addToKey:
						annotKey += c
					else:
						annotValue += c


			currentNodeLen = currentNodeLen.split("[")[0]


	currentNode.name = currentNodeName
	currentNode.length = currentNodeLen

	
	return currentNode


def traverse_depthFirst(currentNode):
	'''
	interates through nodes desendant from currentNode
	'''

	nodeList = [currentNode]

	for child in currentNode.children:
		nodeList += traverse_depthFirst(child)

	#print(nodeList)
	return nodeList



def treeImport_wrap(treefile, inTreeFormat):
	'''
	imports a tree from file and returns parseTree object
	'''
	treeRawString = getTreeRawString(treefile, inTreeFormat)
	openParen_d, closeParn_d, commaParn_d = matchParenCommas(treeRawString)
	root  = buildTree(treeRawString, inTreeFormat, openParen_d, closeParn_d, commaParn_d, 0, len(treeRawString), None)
	t = parseTree(root)

	return t



def makeNewickString(node, treeformat_d):
	'''
	helper for writeNewick which takes current node and recusively builds (chlidren newick)  then appends name as specified in treeformat_d
	'''
	if node.children != []:
		treeString = "("

		for child in node.children:
			childNwk = makeNewickString(child, treeformat_d)
			treeString += childNwk + ","

		treeString = treeString[:-1] + ")"
	else:
		treeString = ""

	treeString += node.name

	if treeformat_d["annotationType"] == "hyphyAnnot":
		treeString += "{"
		for label in node.labels_l:
			treeString += label + ","
		treeString = treeString[:-1] + "}"

	if treeformat_d["branchLen"]:
		treeString += ":" + str(node.length)
		if treeformat_d["annotationType"] == "hyphyAnnot":
			treeString += "[&"
			for annotKey in node.annotKeys:
				annotValue = node.annotKey
				treeString += annotKey + "=" + annotValue + ","
			treeString = treeString[:-1] + "]"



	return treeString


def writeNewick(tree, filename, treeformat):
	'''
	takes a parseTree object and writes in filenames a tree based on specifications of treeformat
	'''

	treeformat_d = parseFormat(treeformat)
	treeString = makeNewickString(tree.root, treeformat_d) + ";"

	directory = os.path.dirname(filename)
	if not os.path.exists(directory):
		os.makedirs(directory)

	file_open = open(filename, "w")
	file_open.write(treeString)
	file_open.close()




	# print("**************Tree: depthfirst*********")
	# for node in traverse_depthFirst(t.root):
	# 	print(node.name, node.length, node.date)#, "date:", node.date, ":", node.l)


