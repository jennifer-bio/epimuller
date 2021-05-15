



class parseTree:
	def __init__(self, root,  name="tree"):
		self.name = name
		self.root = root

class Node:
	def __init__(self, parent, length, name='internalNode'):
		self.parent = parent
		self.name = name
		self.length = length
		self.children = []
		self.annotKeys = []
		
