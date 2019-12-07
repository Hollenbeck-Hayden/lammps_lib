import math_utils as mu


class LammpsDataSet:
	LAMMPS_INT_TYPE = 0
	LAMMPS_FLOAT_TYPE = 1

	TIMESTEP_NAME = 'TimeStep'

	def __init__(self):
		self.data = []
		self.names = {}

	def getNames(self):
		return self.names.keys()
	
	def getType(self, name):
		return self.names[name][0]
	
	def getTypeLen(self, name):
		return self.names[name][1]
	
	def setName(self, name, type, length=1):
		self.names[name] = (type, length)

	def convert(self, value, name):
		if self.getType(name) == LammpsDataSet.LAMMPS_INT_TYPE:
			return int(value)
		elif self.getType(name) == LammpsDataSet.LAMMPS_FLOAT_TYPE:
			return float(value)
		else:
			return value
	
	def checkType(self, name):
		for i in range(len(self.data)):
			if self.getTypeLen(name) == 1:
				self.data[i][name] = self.convert(self.data[i][name], name)
			else:
				self.data[i][name] = [self.convert(value, name) for value in self.data[i][name]]
	
	def checkAll(self):
		for name in self.getNames():
			self.checkType(name)
	
	def addLine(self, line):
		self.data.append(line)

	def select(self, name, index=-1):
		if index < 0:
			return [line[name] for line in self.data]
		elif index < self.getTypeLen(name):
			return [line[name][index] for line in self.data]
		else:
			return []

	def selectID(self, name, id):
		return self.select(name, id-1)
	
	def selectGroup(self, group_name, value, line=0):
		return [self.data[line]['id'][i] for i in range(len(self.data[line][group_name])) if self.data[line][group_name][i] == value]
	
	def selectByID(self, line, name, ids):
		return [self.data[line][name][id-1] for id in ids]
	
	def gather(self, ids, names):
		if type(ids) is not list:
			return [[line[name][ids-1] for name in names] for line in self.data]
		else:
			return [[[line[name][i-1] for name in names] for i in ids] for line in self.data]
	

class LammpsDataFile:
	ATOM_ID = 0
	MOL_ID = 1
	ATOM_TYPE = 2
	ATOM_CHARGE = 3
	ATOM_X = 4
	ATOM_Y = 5
	ATOM_Z = 6

	ATOMS = "atoms"
	BONDS = "bonds"
	ANGLES = "angles"
	DIHEDRALS = "dihedrals"

	def __init__(self):
		self.data = {}
	
	def getNames(self):
		return self.data.keys()
	
	def select(self, name, index=-1, value=-1):
		if value < 0:
			return [line[0] for line in self.data[name]]
		else:
			return [line[0] for line in self.data[name] if line[index] == value]
	
	def get(self, name, i, index):
		return self.data[name][i-1][index]

	def get_atompos(self, atomid):
		return [self.data["atoms"][atomid-1][q] for q in [LammpsDataFile.ATOM_X, LammpsDataFile.ATOM_Y, LammpsDataFile.ATOM_Z]]
	
	def set_atompos(self, atomid, pos):
		self.data["atoms"][atomid-1][LammpsDataFile.ATOM_X] = pos[0]
		self.data["atoms"][atomid-1][LammpsDataFile.ATOM_Y] = pos[1]
		self.data["atoms"][atomid-1][LammpsDataFile.ATOM_Z] = pos[2]

	def displace_atom(self, atomid, dist):
		self.data["atoms"][atomid-1][LammpsDataFile.ATOM_X] += dist[0]
		self.data["atoms"][atomid-1][LammpsDataFile.ATOM_Y] += dist[1]
		self.data["atoms"][atomid-1][LammpsDataFile.ATOM_Z] += dist[2]

	def displace_atoms(self, atomids, dist):
		for atomid in atomids:
			self.displace_atom(atomid, dist)

	def atomids(self):
		return self.select(LammpsDataFile.ATOMS)

	def xcm(self, atoms=[]):
		if len(atoms) == 0:
			atoms = self.atomids()

		total = [0,0,0]
		mass = 0
		for a in atoms:
			v = self.get_atompos(a)
			m = float(self.data["atom types"][self.data["atoms"][a-1][LammpsDataFile.ATOM_TYPE]-1][1])
			mass += m
			for i in range(len(total)):
				total[i] += m * v[i]
		for i in range(len(total)):
			total[i] /= mass
		return total
	
	def get_atom_ids(self):
		return [line[LammpsDataFile.ATOM_ID] for line in self.data["atoms"]]
	
	def add(self, name, line):
		if name not in self.getNames():
			self.data[name] = []
		line.insert(0, len(self.data[name])+1)
		self.data[name].append(line)
		return line[0]

	def addType(self, name, line):
		return self.add(LammpsDataFile.getTypeName(name), line)
	
	def organize(self):
		for name in self.getNames():
			self.data[name] = sorted(self.data[name], key=lambda x: x[0])
	
	@staticmethod
	def getTypeName(name):
		return name[:-1] + " types"

	def addAtom(self, molid, atomtype, charge, position):
		return self.add(LammpsDataFile.ATOMS, [molid, atomtype, charge, position[0], position[1], position[2]])

	def addBond(self, bondtype, atom1, atom2):
		return self.add(LammpsDataFile.BONDS, [bondtype, atom1, atom2])
	
	def addAngle(self, angletype, atom1, atom2, atom3):
		return self.add(LammpsDataFile.ANGLES, [angletype, atom1, atom2, atom3])
	
	def addDihedral(self, dihedraltype, atom1, atom2, atom3, atom4):
		return self.add(LammpsDataFile.DIHEDRALS, [dihedraltype, atom1, atom2, atom3, atom4])

	# CURRENTLY DOESN'T SUPPORT BOX SIZE
	def writeToFile(self, filename):
		outfile = open(filename, "w")
		names = ["atom", "bond", "angle", "dihedral", "improper"]
		# Write header
		outfile.write("LAMMPS Description\n\n")
		for n in names:
			name = n + "s"
			if name not in self.getNames():
				continue
			outfile.write("\t" + str(len(self.data[name])) + " " + name + "\n")

		for n in names:
			name = n + " types"
			if name not in self.getNames():
				continue
			outfile.write("\t" + str(len(self.data[name])) + " " + name + "\n")

		outfile.write("\t-100 100 xlo xhi\n")
		outfile.write("\t-100 100 ylo yhi\n")
		outfile.write("\t-100 100 zlo zhi\n")

		# Write data
		outfile.write("Masses\n\n")
		for line in self.data["atom types"]:
			for item in line:
				outfile.write(str(item))
				outfile.write(" ")
			outfile.write("\n")

		for n in names:
			name = n + "s"
			if name not in self.getNames():
				continue
			outname = name[0].upper() + name[1:]
			outfile.write("\n" + outname + "\n\n")
			for line in self.data[name]:
				for item in line:
					outfile.write(str(item))
					outfile.write(" ")
				outfile.write("\n")

		for n in names[1:]:
			name = n + " types"
			if name not in self.getNames():
				continue
			outname = n[0].upper() + n[1:] + " Coeffs"
			outfile.write("\n" + outname + "\n\n")
			for line in self.data[name]:
				for item in line:
					outfile.write(str(item))
					outfile.write(" ")
				outfile.write("\n")

		outfile.close()


def readFixFile(filename):
	infile = open(filename, "r")
	dataset = LammpsDataSet()
	name_order = []
	for line in infile:
		# Read names
		if line.startswith("# TimeStep"):
			scalars = []
			items = line.strip().split()
			for i in range(1, len(items)):
				name = items[i]

				# Determine type length
				typelen = 1
				if name.endswith("]"):
					name = items[i].split("[")[0]
					typelen = int(items[i].split("[")[1][:-1])
					if name in scalars:
						typelen += 1
				else:
					scalars.append(name)

				# Determine type
				type = LammpsDataSet.LAMMPS_FLOAT_TYPE
				if name == LammpsDataSet.TIMESTEP_NAME:
					type = LammpsDataSet.LAMMPS_INT_TYPE

				dataset.setName(name, type, typelen)
				if name not in name_order: name_order.append(name)
					

		# Ignore comments
		elif line.startswith("#"):
			continue

		# Read data
		else:
			data = {}
			items = line.strip().split()
			index = 0
			for i in range(len(name_order)):
				name = name_order[i]
				temp = []
				if dataset.getTypeLen(name) > 1:
					for j in range(dataset.getTypeLen(name)):
						temp.append(dataset.convert(items[index+j], name))
				else:
					temp = dataset.convert(items[index], name)
				data[name] = temp
				index += dataset.getTypeLen(name)
			dataset.addLine(data)
	infile.close()
	return dataset


def readDumpFile(filename):
	infile = open(filename, "r")
	dataset = LammpsDataSet()
	NONE = 0
	READ_TIMESTEP = 1
	READ_NUM_ATOMS = 2
	READ_BOX_BOUNDS = 3
	READ_ATOMS = 4
	state = NONE
	data = {}
	name_order = []
	for line in infile:
		if line.startswith("ITEM: TIMESTEP"):
			state = READ_TIMESTEP
			if len(data.keys()) > 0:
				dataset.addLine(data)
			data = {}
			name_order = []

		elif line.startswith("ITEM: NUMBER OF ATOMS"):
			state = READ_NUM_ATOMS

		elif line.startswith("ITEM: BOX BOUNDS"):
			state = READ_BOX_BOUNDS

		elif line.startswith("ITEM: ATOMS"):
			state = READ_ATOMS
			items = line.strip().split()
			for name in items[2:]:
				type = LammpsDataSet.LAMMPS_FLOAT_TYPE
				if name == "id":
					type = LammpsDataSet.LAMMPS_INT_TYPE
				elif name == "type":
					type = LammpsDataSet.LAMMPS_INT_TYPE
				elif name == "mol":
					type = LammpsDataSet.LAMMPS_INT_TYPE
				dataset.setName(name, type, data["NumAtoms"])
				if name not in name_order: name_order.append(name)
				data[name] = [0 for i in range(data["NumAtoms"])]

		elif state == READ_TIMESTEP:
			data[LammpsDataSet.TIMESTEP_NAME] = int(line.strip())
			dataset.setName(LammpsDataSet.TIMESTEP_NAME, LammpsDataSet.LAMMPS_INT_TYPE, 1)
			state = NONE

		elif state == READ_NUM_ATOMS:
			data["NumAtoms"] = int(line.strip())
			dataset.setName("NumAtoms", LammpsDataSet.LAMMPS_INT_TYPE, 1)
			state = NONE

		elif state == READ_BOX_BOUNDS:
			pass

		elif state == READ_ATOMS:
			items = line.strip().split()
			index = int(items[0])-1
			for i in range(len(items)):
				data[name_order[i]][index] = dataset.convert(items[i], name)

		else:
			pass

	infile.close()
	return dataset

def readDataFile(filename):
	infile = open(filename, "r")
	dataset = LammpsDataFile()
	state = ""
	for line in infile:
		# remove comments
		if "#" in line:
			index = line.find("#")
			line = line[:index]

		line = line.strip()

		if len(line) == 0:
			continue
		elif line == "LAMMPS Description":
			state = "desc"
			continue
		elif line == "Masses":
			state = "atom types"
			continue

		items = line.split()
		if not items[0].isdigit():
			state = line.lower()
			continue
		
		if state == "desc":
			if len(items) == 2:
				name = items[1].lower()
			elif len(items) == 3:
				name = (items[1] + " " + items[2]).lower()
			else:
				continue

			count = int(items[0])
			dataset.data[name] = [0 for i in range(count)]
		elif state == "atoms":
			dataset.data[state][int(items[0])-1] = [int(items[0]), int(items[1]), int(items[2]), float(items[3]),
						float(items[4]), float(items[5]), float(items[6])]
		elif state.endswith("types"):
			i = int(items[0])-1
			dataset.data[state][i] = [float(x) for x in items[1:]]
			dataset.data[state][i].insert(0, i+1)
		else:
			i = int(items[0])-1
			dataset.data[state][i] = [int(x) for x in items]
	infile.close()
	return dataset



# Computes the average of a list
def averageList(l):
	total = 0
	for i in l:
		total += i
	return total / len(l)

# Combines datasets if they share timesteps
def combineDatasets(d1, d2):
	dataset = LammpsDataSet()

	for name in d1.getNames(): dataset.setName(name, d1.getType(name), d1.getTypeLen(name))
	for name in d2.getNames(): dataset.setName(name, d2.getType(name), d2.getTypeLen(name))

	i1 = 0
	i2 = 0
	while i1 < len(d1.data) and i2 < len(d2.data):
		time1 = d1.data[i1][LammpsDataSet.TIMESTEP_NAME]
		time2 = d2.data[i2][LammpsDataSet.TIMESTEP_NAME]
		if time1 == time2:
			line = {}
			for name in d1.getNames(): line[name] = d1.data[i1][name]
			for name in d2.getNames(): line[name] = d2.data[i2][name]
			dataset.addLine(line)
			i1 += 1
		elif time1 < time2:
			i1 += 1
		else:
			i2 += 1
	
	return dataset


def distFilename(dist):
	return "{0:.1f}".format(dist)
