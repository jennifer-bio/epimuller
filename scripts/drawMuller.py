import fnmatch
import os
import argparse
import sys
import math
import svgwrite #conda install svgwrite
import random
import cairo
import cairosvg


# WIDTH = 1500
# HEIGHT = 1000
# LEGENDWIDTH = 220
# MARGIN = 60
# LABELSHIFT = 15
# FONTSIZE = 26

class Clade:
	def __init__(self, name, parent_name):
		self.name = name
		self.parent_name = parent_name
		self.cladeSnapshot_time_d = {} 
		self.y1_d = {}
		self.y2_d = {}

class Snapshot:
	def __init__(self, time, date):
		self.time = time
		self.date = date
		self.label = "t_" + str(time)
		self.cladeSnapshot_clade_d = {}
		self.sumAll = 0 


class CladeSnapshot:
	def __init__(self, clade, snapshot, abundance):
		self.clade = clade
		self.snapshot = snapshot
		self.abundance = abundance
		self.sumDescendant = 0 
		#self.numChildren =  len(self.clade.children)
		#self.numChildren = 0
		#to be added later sumDescendant , numChildren

	def sumUpDescendants(self):
		summed = 0
		for child in self.clade.children:
			child_tSnap = child.cladeSnapshot_time_d[self.snapshot.time]
			if child.children_count > 0:
				summed += child_tSnap.sumUpDescendants()
			summed += child_tSnap.abundance
		return summed


def makeColor():
	# random_number = random. randint(0,16777215)
	# color = "#" + str(hex(random_number))[2:]
	color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
	return color

def timeToX(time, scaleTime, minTime):
	return(scaleTime*(int(time) - int(minTime)) + MARGIN)



def textwidth(text, fontsize):
	#function copied from http://blog.mathieu-leplatre.info/text-extents-with-python-cairo.html
	try:
		import cairo
	except Exception as e:
		return len(str) * fontsize
	surface = cairo.SVGSurface('undefined.svg', 1280, 200)
	cr = cairo.Context(surface)
	cr.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
	cr.set_font_size(fontsize)
	xbearing, ybearing, width, height, xadvance, yadvance = cr.text_extents(text)
	return width


def textheight(text, fontsize):
	#function based on textwidth copied from http://blog.mathieu-leplatre.info/text-extents-with-python-cairo.html
	try:
		import cairo
	except Exception as e:
		return len(str) * fontsize
	surface = cairo.SVGSurface('undefined.svg', 1280, 200)
	cr = cairo.Context(surface)
	cr.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
	cr.set_font_size(fontsize)
	xbearing, ybearing, width, height, xadvance, yadvance = cr.text_extents(text)
	return height

# def makeNoScaleCord(abundances_d, childParent_d, timeToDate_d):
# 	'''
# 	returns list of shape_l where shape_l is a list cord_l, which is cordinates for a clade, fomated as [(x, y), (x, y), ...] relative to top left
# 	'''
# 	shape_l	= []
# 	return shape_l	

def drawWrapper(outFolder, outPrefix, root_clades_l, scaleTime, times_l, maxY, minTime, labelPosition, xlabel, timeToDate_d):

	#Draw background

	outFile =  os.path.join(outFolder, outPrefix + ".svg")
	outFilePDF =  os.path.join(outFolder, outPrefix + ".pdf")

	img = svgwrite.Drawing(filename = outFile, size = (str(WIDTH)+"px", str(HEIGHT)+"px"))
	img.add(img.polyline(points = [(0,0), (0, HEIGHT), (WIDTH, HEIGHT), (WIDTH, 0)], stroke='black', fill = 'white'))
	#img.add(img.text(text = 'Legend', insert = (WIDTH-LEGENDWIDTH, LEGENDWIDTH), font_size=24))
	#rightWidth = (WIDTH-(MARGIN + LEGENDWIDTH + 1/scaleTime)) #WIDTH - (MARGIN + LEGENDWIDTH) 
	rightWidth = timeToX(list(times_l)[-1], scaleTime, minTime)
	img.add(img.polyline(points = [(MARGIN,MARGIN), (MARGIN, HEIGHT-MARGIN), (rightWidth, HEIGHT-MARGIN), (rightWidth, MARGIN), (MARGIN,MARGIN)], stroke='black', stroke_width=5, fill="white"))


	#draw clades
	img, x_labelCord_l, y_labelCord_l, label_l = extractCord_draw(root_clades_l, img, scaleTime, [], [], [], times_l, minTime, labelPosition )


	tHeight = textheight("TESTLABEL", FONTSIZE)
	#write clade labels
	if labelPosition == "Right":
		xpos = timeToX(list(times_l)[-1], scaleTime, minTime)+LABELSHIFT

		zipped = zip(label_l, y_labelCord_l)
		sort_zip = list(sorted(zipped, key = lambda x: x[1]))

		tWidth = textwidth(sort_zip[0][0], FONTSIZE)
		fontsize = FONTSIZE 
		while tWidth > LEGENDWIDTH:
			fontsize += (-1)
			tWidth = textwidth(sort_zip[0][0], fontsize)

		ypos = sort_zip[0][1]
		img.add(img.text(text = sort_zip[0][0], insert = (xpos, ypos), font_size=fontsize))
		for i in range(1, len(sort_zip)):
			bottomPrevious = ypos + LABELSHIFT + tHeight
			ypos = sort_zip[i][1]
			if ypos < bottomPrevious:
				ypos = bottomPrevious

			tWidth = textwidth(sort_zip[i][0], FONTSIZE)
			fontsize = FONTSIZE 
			while tWidth > (LEGENDWIDTH+MARGIN):
				fontsize += (-1)
				tWidth = textwidth(sort_zip[i][0], fontsize)
			
			img.add(img.text(text = sort_zip[i][0], insert = (xpos, ypos), font_size=fontsize))
	else: #if args.labelPosition == "Start":
		for i in range(len(label_l)):
			img.add(img.text(text = label_l[i], insert = (x_labelCord_l[i], y_labelCord_l[i]), font_size=FONTSIZE))


	#wirte y axis labels
	topPlot = MARGIN
	bottomPlot = HEIGHT-MARGIN
	totalHeight = bottomPlot-topPlot
	wirteEvery = int(maxY/10)
	maxLabel = str(int(maxY))+ " "
	tWidth = textwidth(maxLabel, FONTSIZE)
	fontsize = FONTSIZE 
	while tWidth > MARGIN:
		fontsize = fontsize - 1 
		tWidth = textwidth(maxLabel, fontsize)
	for i in range(int(maxY)):
		if i % wirteEvery == 0:
			xpos = MARGIN - (textwidth(str(i), fontsize) + LABELSHIFT)
			ypos = MARGIN + (maxY - i)*totalHeight/maxY
			img.add(img.text(text = str(i), insert = (xpos,  ypos + (textheight(str(i), FONTSIZE))/2), font_size=fontsize))
			img.add(img.line(start = (MARGIN-LABELSHIFT,  ypos), end = (MARGIN+LABELSHIFT,  ypos), stroke_width=5, stroke = "black"))




	#write x axis labels
	if xlabel == "date":
		tWidth = textwidth("2021-05-03", FONTSIZE)
	else:
		tWidth = textwidth("200", FONTSIZE)
	numLabels = rightWidth/tWidth 
	wirteEvery = math.ceil(len(times_l)/numLabels)
	for time in times_l:
		if int(time) % wirteEvery == 0:
			#if int(time) % int(args.XLABFREQ) == 0:
			if xlabel == "time":
				label = str(time)
			else:
				label = timeToDate_d[time]
			img.add(img.text(text = label, insert = (timeToX(time, scaleTime, minTime),  HEIGHT-MARGIN+(2*LABELSHIFT)), font_size=FONTSIZE))
			img.add(img.line(start = (timeToX(time, scaleTime, minTime),  HEIGHT-MARGIN+LABELSHIFT), end = (timeToX(time, scaleTime, minTime),  HEIGHT-(MARGIN+LABELSHIFT)), stroke_width=5, stroke = "black"))
	
	# write title
	fontsize = FONTSIZE+6
	tHeight = textheight(outPrefix, fontsize)
	while tHeight/2 > MARGIN:
		fontsize = fontsize - 1 
		tHeight = textwidth(outPrefix, fontsize)

	img.add(img.text(text = outPrefix.replace("_", " "), insert = (MARGIN + LABELSHIFT, MARGIN/2), font_size=fontsize))

	img.save()

	cairosvg.svg2pdf(file_obj=open(outFile, "rb"), write_to=outFilePDF)

	




def extractCord_draw(clades_l, img, scaleTime, x_labelCord_l, y_labelCord_l, label_l, times_l, minTime, labelPosition ):
	#for snap in allTimes_l:
	clade_cord_d = {} #key clade name, value is list of coordinate tuples [... (x2, y2t2), (x1, y2t1), (x1, y1t1), (x2, y1t2) ...]
	#clade_col_d = {} #key clade name, value is color
	for clade in clades_l:
		coordinate_l = []
		startAbundance_time = 'NA'
		lastAbundance_time = 0
		maxAbundance_time = 'NA'
		maxAbundance_value = -1
		cladeDrawn = False
		for time in times_l:
			#x = int(time)*scaleTime
			if time in clade.y1_d:
				y1 = clade.y1_d[time]
				y2 = clade.y2_d[time]
				coordinate_l = [(timeToX(time, scaleTime, minTime), y1)] + coordinate_l +  [(timeToX(time, scaleTime, minTime), y2)]
				abundance = clade.cladeSnapshot_time_d[time].abundance/(1.0*clade.cladeSnapshot_time_d[time].snapshot.sumAll)

				if abundance > 0:
					cladeDrawn = True
					if startAbundance_time == 'NA':
						startAbundance_time = time
						startAbundance_y = (y2-y1)//2 + y1
					lastAbundance_time = time
					lastAbundance_y = y1
					drawAbundance = y2-y1
					if maxAbundance_value < drawAbundance: #abundance:
						maxAbundance_value = drawAbundance
						maxAbundance_time = time
						maxAbundance_y = drawAbundance//2 + y1

		clade_cord_d[clade.name] = coordinate_l

		if cladeDrawn:
			if labelPosition == "Right":
				x_cord = timeToX(list(times_l)[-1], scaleTime, minTime)+LABELSHIFT
				y1_cord = y1+LABELSHIFT
			elif labelPosition == "Max":
				x_cord = timeToX(maxAbundance_time, scaleTime, minTime) - textwidth(clade.name, FONTSIZE)//2
				if x_cord < 0:
					x_cord = 0
				y1_cord = maxAbundance_y
			elif labelPosition == "Start":
				x_cord = timeToX(startAbundance_time, scaleTime, minTime)-(textwidth(clade.name, FONTSIZE))
				y1_cord = startAbundance_y
				if x_cord < 0:
					x_cord = 0
			elif labelPosition == "End":
				x_cord = timeToX(lastAbundance_time, scaleTime, minTime)+LABELSHIFT
				if x_cord < 0:
					x_cord = 0
				y1_cord = lastAbundance_y
			else:
				print("do not know where to put clade label")
				sys.exit(1)
			x_labelCord_l.append(x_cord)
			y_labelCord_l.append(y1_cord)


			label_l.append(clade.name)

		
		img.add(img.polyline(points = coordinate_l, stroke='black', stroke_width=0, fill=clade.color))

		img, x_labelCord_l, y_labelCord_l, label_l = extractCord_draw(clade.children, img, scaleTime, x_labelCord_l, y_labelCord_l, label_l, times_l, minTime, labelPosition )


	return(img, x_labelCord_l, y_labelCord_l, label_l)

def removeSmallClades(abundances_d, heiarchy_d, minCount):
	'''
	
	removes clade from abudances and heiarchiy that have sum less than minCount
	heiarchy_d: key:child clade; value:parent clade
	abundances_d: key: week; value: dict of key:clade; value: count
	'''

	# remove clades with less than minCount 

	avalibleClades_d = {}
	avalibleClades_d['anc'] = 0
	avalibleClades_d['NA'] = 0

	avalibleClades_s = set()
	removedClades_s = set()

	for week in abundances_d:
		for clade in abundances_d[week]:
			if clade not in avalibleClades_d:
				avalibleClades_d[clade] = 0
			avalibleClades_d[clade] += int(abundances_d[week][clade])

	old_heiarchy_d = heiarchy_d.copy()
	for clade in old_heiarchy_d: #old_heiarchy_d.keys():
		if clade not in avalibleClades_d:
			heiarchy_d.pop(clade)
			removedClades_s.add(clade)
		elif avalibleClades_d[clade] < int(minCount) and clade != 'NA' and clade != 'anc':
			heiarchy_d.pop(clade)
			removedClades_s.add(clade)
		else:
			avalibleClades_s.add(clade)

	for child in heiarchy_d:
		parent = old_heiarchy_d[child]
		while parent not in heiarchy_d and parent != 'NA':
			newParent = old_heiarchy_d[parent] #heiarchy_d[parent]
			parent = newParent
		heiarchy_d[child] = parent

	# remove from abundances if remved from hieracrhy

	old_abundances_d = abundances_d.copy()
	for week in list(old_abundances_d):
		for clade in list(old_abundances_d[week]):
			if clade in removedClades_s:
				abundances_d[week].pop(clade)

	return(abundances_d, heiarchy_d)


def defineChildBoundries(time, scaleFactor, parentCladeSnap, y1_parent, y2_parent):
	'''
	time is the number defining the time of interest
	parentCladeSnap is pointer to snapshotClade object which has boundries of"
	y1_parent is top of clade boundries, y1_parent is bottom of clade boundries

	updates 
	'''
	numSections = len(parentCladeSnap.clade.children) + 1
	#descendantSpace = (parentCladeSnap.sumDescendant - parentCladeSnap.abundance)
	#parentEdgeHeight = scaleFactor*(parentCladeSnap.abundance - descendantSpace)/numSections	
	parentEdgeHeight = scaleFactor*(parentCladeSnap.abundance)/numSections

	acountedHeight = y1_parent
	for childClade in parentCladeSnap.clade.children:
		childCladeSnap = childClade.cladeSnapshot_time_d[time]
		y1 =  acountedHeight + parentEdgeHeight
		y2 = y1 + (scaleFactor*childCladeSnap.sumDescendant)
		acountedHeight += parentEdgeHeight + (scaleFactor*childCladeSnap.sumDescendant)
		childClade.y1_d[time] = y1
		childClade.y2_d[time] = y2
		defineChildBoundries(time, scaleFactor, childCladeSnap, y1, y2)


def main():
##########################  parse user arguments
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)



	parser.add_argument('-p', '--parentHierarchy_name', required=True, type=str, help="csv output from mutationLinages_report.py with child parent col")

	parser.add_argument('-a', '--abundance_name', required=True, type=str, help="csv output from mutationLinages_report.py with abundances of clades")

	parser.add_argument('-c', '--cases_name', required=False, type=str, help="file with cases - formated with 'date' in ISO format and 'confirmed_rolling' cases, in tsv format")

	parser.add_argument('-o', '--outFolder', required=True, type=str, help="csv output from mutationLinages_report.py with child parent col") 


	parser.add_argument('-mt', '--MINTIME', required=False, type=str, default="30", help="minimum time point to start plotting")
	parser.add_argument('-min', '--MINTOTALCOUNT', required=False, type=str, default="50", help="minimum total count for group to be included")

	parser.add_argument('-l', '--xlabel', required=False, type=str, choices = ["date", "time"], default="date", help="Format of x axis label: ISO date format or timepoints from start")
	parser.add_argument('-lp', '--labelPosition', required=False, type=str, default="Right", choices = ["Right", "Max", "Start", "End"], help="choose position of clade labels")

	#parser.add_argument("--tmrca", action="store_true", help="draw point at tmrca of clade if flag is used")


	drawing_group_page  = parser.add_argument_group('Options for page setup')

	drawing_group_page.add_argument('--WIDTH', required=False, type=str, default="1500", help="WIDTH of page (px)")
	drawing_group_page.add_argument('--HEIGHT', required=False, type=str, default="1000", help="HEIGHT of page (px)")
	drawing_group_page.add_argument('--LEGENDWIDTH', required=False, type=str, default="220", help="LEGENDWIDTH to the right of plotting area (px)")
	drawing_group_page.add_argument('--LABELSHIFT', required=False, type=str, default="15", help="nudge label over by LABELSHIFT (px)")
	drawing_group_page.add_argument('--MARGIN', required=False, type=str, default="60", help="MARGIN around all sides of plotting area (px)")
	drawing_group_page.add_argument('--FONTSIZE', required=False, type=str, default="26")
	

	args = parser.parse_args()


########################## set up global


	if not os.path.exists(args.outFolder):
		os.makedirs(args.outFolder)

	global WIDTH
	WIDTH = float(args.WIDTH)
	global HEIGHT
	HEIGHT = float(args.HEIGHT)
	global LEGENDWIDTH
	LEGENDWIDTH = float(args.LEGENDWIDTH)
	global MARGIN
	MARGIN = float(args.MARGIN)
	global FONTSIZE
	FONTSIZE = float(args.FONTSIZE)
	global LABELSHIFT
	LABELSHIFT = float(args.LABELSHIFT)

	if args.labelPosition != "Right":
		LEGENDWIDTH = 0
	

########################## read in files

	abundances_file = open(args.abundance_name, "r")

	abundances_d = {}
	#abundance_total_d {}
	timeToDate_d = {}
	times_index = 'NA'
	for line in abundances_file:
		line_l = line.strip().split(",")
		if line_l[0] == "names":
			times_index = line_l.index("times")
			a_index = line_l.index("abundances")
			date_index = line_l.index("date")
		elif times_index == 'NA':
			print("First line of abundances_file must start with 'names' col and contain 'times', 'abundances', and 'date' cols, with no spaces between commas\n")
			sys.exit(1)
		else:
			time = line_l[times_index]
			if int(time) >= int(args.MINTIME):
				clade = line_l[0]
				abundance = line_l[a_index]
				if time not in abundances_d:
					abundances_d[time] = {}
				#if clade not in abundance_total_d:
				#	abundance_total_d[clade] = 0
				#abudance_total_d[clade] += abundance
				abundances_d[time][clade] = abundance
				if time not in timeToDate_d:
					timeToDate_d[time] = line_l[date_index]

	abundances_file.close()	


	hierarchy_file = open(args.parentHierarchy_name, "r")	

	childParent_d = {}
	cladeColor_d = {}
	hasHeader = False
	hasColor = False
	for line in hierarchy_file:
		line_l = line.strip().split(",")
		if line_l[0] == "names" and line_l[1] == "parents":
			hasHeader = True
			if len(line_l) > 2:
				if line_l[2] == "color":
					hasColor = True
		elif not hasHeader:
			print("First line of parentHierarchy_name must have 'names' as first col and 'parents' as second col, sperated with commas and no spaces")
			sys.exit(1)
		else:
			childParent_d[line_l[0]] = line_l[1]
			if hasColor:
				if line_l[2] != "":
					cladeColor_d[line_l[0]] = line_l[2]
				else:
					cladeColor_d[line_l[0]] = makeColor()
			else:
				cladeColor_d[line_l[0]] = makeColor()

	hierarchy_file.close()

	abundances_d, heiarchy_d = removeSmallClades(abundances_d, childParent_d, args.MINTOTALCOUNT)


########################## parse clades
	times_l = abundances_d.keys()
	#timeLabs_s = set()
	clades_l = []
	root_clades_l = []

	for clade in childParent_d:
		clade = Clade(clade, childParent_d[clade])
		setattr(clade, "color", cladeColor_d[clade.name])
		#TODO make this determatistic for most diverse colors and add in more color options
		clades_l.append(clade)
		if clade.parent_name == "NA":
			root_clades_l.append(clade)


	#add pointers to children clade objects
	for parent in clades_l:
		setattr(parent, "children", [])
		for child in clades_l:
			if child.parent_name == parent.name:
				parent.children = parent.children + [child]
		setattr(parent, "children_count", len(parent.children))

	#add pointers to parent clade objects
	for child in clades_l:
		setattr(child, "parent", Clade("NA", "NA"))
		for parent in clades_l:
			if parent.name == child.parent_name:
				setattr(child, "parent", parent)
			

	#record abundances in objects
	allTimes_l = []
	for time in abundances_d:
		if time >= args.MINTIME:
			t = Snapshot(time, timeToDate_d[time])
			for clade in clades_l:
				if clade.name in abundances_d[time]:
					abundance = int(abundances_d[time][clade.name])
				else:
					abundance = 0
				clade_oneTime = CladeSnapshot(clade, t, abundance)
				t.cladeSnapshot_clade_d[clade.name] = clade_oneTime
				t.sumAll += abundance
				clade.cladeSnapshot_time_d[time] = clade_oneTime
			allTimes_l.append(t)

	#calculate Descendants and children
	for snap in allTimes_l:
		for cladeSnap in snap.cladeSnapshot_clade_d.values():
			cladeSnap.sumDescendant = cladeSnap.abundance + cladeSnap.sumUpDescendants()

# ################################ determain plotting values

	
	topPlot = MARGIN
	bottomPlot = HEIGHT-MARGIN

	for snap in allTimes_l:
		totalHeight = (bottomPlot-topPlot)
		scaleFactor = totalHeight/snap.sumAll

		descendantSpace = 0
		for clade in root_clades_l:
			descendantSpace += clade.cladeSnapshot_time_d[snap.time].sumDescendant

		numSections = len(root_clades_l) + 1
		parentEdgeHeight = scaleFactor*(snap.sumAll - descendantSpace)/numSections
		acountedHeight = topPlot


		for clade in root_clades_l:
			cladeSnap = clade.cladeSnapshot_time_d[snap.time]
			y1 = acountedHeight + parentEdgeHeight
			y2 = y1 + (scaleFactor*cladeSnap.sumDescendant)
			acountedHeight += (scaleFactor*cladeSnap.sumDescendant) + parentEdgeHeight
			clade.y1_d[snap.time] = y1
			clade.y2_d[snap.time] = y2
			defineChildBoundries(snap.time, scaleFactor, cladeSnap, y1, y2)


	scaleTime = (WIDTH-(MARGIN+LEGENDWIDTH))/len(times_l)
	drawWrapper(args.outFolder, "relative_abundance", root_clades_l, scaleTime, times_l, 100, args.MINTIME, args.labelPosition, args.xlabel, timeToDate_d)



# ########################## make fig with number of samples scaling 

	
	topPlot = MARGIN
	bottomPlot = HEIGHT-MARGIN

	maxCount = 0
	for snap in allTimes_l:
		if snap.sumAll > maxCount:
			maxCount = snap.sumAll


	for snap in allTimes_l:
		#adjust relative to scalling
		lessThanMax = maxCount - snap.sumAll
		lessThanMax_ratio = snap.sumAll/maxCount
		totalHeight = (bottomPlot-topPlot)*lessThanMax_ratio
		
		scaleFactor = totalHeight/snap.sumAll 

		descendantSpace = 0
		for clade in root_clades_l:
			descendantSpace += clade.cladeSnapshot_time_d[snap.time].sumDescendant

		numSections = len(root_clades_l) + 1
		parentEdgeHeight = scaleFactor*(snap.sumAll - descendantSpace)/numSections
		acountedHeight = topPlot + scaleFactor*lessThanMax

		for clade in root_clades_l:

			cladeSnap = clade.cladeSnapshot_time_d[snap.time]
			y1 = acountedHeight + parentEdgeHeight
			y2 = y1 + (scaleFactor*cladeSnap.sumDescendant)
			acountedHeight += (scaleFactor*cladeSnap.sumDescendant) + parentEdgeHeight
			clade.y1_d[snap.time] = y1
			clade.y2_d[snap.time] = y2
			defineChildBoundries(snap.time, scaleFactor, cladeSnap, y1, y2)


	scaleTime = (WIDTH-(MARGIN+LEGENDWIDTH))/len(times_l)
	drawWrapper(args.outFolder, "sequence_scaled_lineages", root_clades_l, scaleTime, times_l, maxCount,   args.MINTIME, args.labelPosition, args.xlabel, timeToDate_d)



########################## make fig with cases scaling
	if args.cases_name is not None:

		if "strain" in line_l and "date" in line_l:
			strain_index = line_l.index("strain")
			date_index = line_l.index("date")

		date_index = "na"
		case_index = "na"
		dateToCase_d = {}
		cases_file = open(args.cases_name, "r")
		for line in cases_file:
			line_l = line.strip().split("	")
			if date_index == "na":
				if "confirmed_rolling" in line_l and "date" in line_l:
					date_index = line_l.index("date")
					case_index = line_l.index("confirmed_rolling")
				else:
					print("file with cases - formated with 'date' in ISO format and 'confirmed_rolling' cases, in tsv format")
					sys.exit(1)
			else:
				dateToCase_d[line_l[date_index].replace('"', "")] = float(line_l[case_index].replace('"', ""))


		topPlot = MARGIN
		bottomPlot = HEIGHT-MARGIN

		maxCount = 0
		for snap in allTimes_l:
			rollCase = dateToCase_d[timeToDate_d[snap.time]]
			if rollCase > maxCount:
				maxCount = rollCase

		for snap in allTimes_l:

			rollCase = dateToCase_d[timeToDate_d[snap.time]]
			cases_lessThanMax = maxCount - rollCase
			psudoAbundance_lessThanMax = (cases_lessThanMax/rollCase)*snap.sumAll
			lessThanMax_ratio = rollCase/maxCount
			totalHeight = (bottomPlot-topPlot)*lessThanMax_ratio

			
			scaleFactor = totalHeight/snap.sumAll 


			descendantSpace = 0
			for clade in root_clades_l:
				descendantSpace += clade.cladeSnapshot_time_d[snap.time].sumDescendant

			numSections = len(root_clades_l) + 1
			parentEdgeHeight = scaleFactor*(snap.sumAll - descendantSpace)/numSections
			acountedHeight = topPlot + psudoAbundance_lessThanMax*scaleFactor


			for clade in root_clades_l:

				cladeSnap = clade.cladeSnapshot_time_d[snap.time]
				y1 = acountedHeight + parentEdgeHeight
				y2 = y1 + (scaleFactor*cladeSnap.sumDescendant)
				acountedHeight += (scaleFactor*cladeSnap.sumDescendant) + parentEdgeHeight
				clade.y1_d[snap.time] = y1
				clade.y2_d[snap.time] = y2
				defineChildBoundries(snap.time, scaleFactor, cladeSnap, y1, y2)


		scaleTime = (WIDTH-(MARGIN+LEGENDWIDTH))/len(times_l)
		drawWrapper(args.outFolder, "case_scaled_lineages", root_clades_l, scaleTime, times_l, maxCount, args.MINTIME, args.labelPosition, args.xlabel, timeToDate_d)
	else:
		print("No case data supplied - skipping case scaled plot")



if __name__ == "__main__":
	main()