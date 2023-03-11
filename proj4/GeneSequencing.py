#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random


# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):

		self.banded = banded
		self.MaxCharactersToAlign = align_length

		#Debugging
		#seq1 = 'exponential'
		#seq2 = 'polynomial'

		delin = 5
		sub = 1
		found = -3

		if (banded and (align_length > len(seq1) or align_length > len(seq2))):
			return self.align_banded(seq1, seq2, align_length)

		#print('starting align ************************************')
		#print(seq1)
		#print(seq2)
		left = 1
		up = 2
		diagonal = 3

		costArray = []
		prevArray = []
		array = []
		for i in range(len(seq1) + 1):
			array.append(0)
		for j in range(len(seq2) + 1):
			costArray.append(array.copy())
			prevArray.append(array.copy())
		print(costArray)

		#populate array
		print('populating-------------------------------')
		for j in range(len(seq2) + 1):
			costArray[j][0] = j*delin
			if (j != 0):
				prevArray[j][0] = up
		for i in range(len(seq1) + 1):
			costArray[0][i] = i*delin
			if (i != 0):
				prevArray[0][i] = left
		print('next')
		print(costArray)
		for j in range(1, len(costArray)):
			for i in range(1, len(costArray[0])):
				subCost = costArray[j - 1][i - 1] + sub
				if (seq1[i-1] == seq2[j-1]):
					subCost = subCost - sub + found
				costArray[j][i] = min(subCost, costArray[j-1][i] + delin, costArray[j][i-1] + delin)

				#go in order left, up, diagonal to determine what previous is assigned
				if (costArray[j][i] == costArray[j][i-1] + delin):
					prevArray[j][i] = left
				elif (costArray[j][i] == costArray[j-1][i] + delin):
					prevArray[j][i] = up
				else:
					prevArray[j][i] = diagonal

		#debug
		#print('results---------------------------------')
		#print('cost:')
		#print(costArray)
		#print('prev')
		#print(prevArray)

		#getting the alignment
		i = len(prevArray[0]) - 1
		j = len(prevArray) - 1
		alignSeq1 = ''
		alignSeq2 = ''
		while (prevArray[j][i] != 0):
			if (prevArray[j][i] == left):
				alignSeq1 = seq1[i-1] + alignSeq1
				alignSeq2 = '-' + alignSeq2
				i -= 1
			elif (prevArray[j][i] == up):
				alignSeq1 = '-' + alignSeq1
				alignSeq2 = seq2[j-1] + alignSeq2
				j -= 1
			elif (prevArray[j][i] == diagonal):
				alignSeq1 = seq1[i-1] + alignSeq1
				alignSeq2 = seq2[j-1] + alignSeq2
				i -= 1
				j -= 1
			
			#trim so that only the first 100 characters are visible
			if (len(alignSeq1) > 100):
				alignSeq1 = alignSeq1[0:len(alignSeq1) - 1]
			if (len(alignSeq2) > 100):
				alignSeq2 = alignSeq2[0:len(alignSeq2) - 1]

		print('alignment------------')
		print(alignSeq1)
		print(alignSeq2)
		print(len(alignSeq1))
		print(len(alignSeq2))
		print('cost-----------------')
		print(costArray[len(costArray) - 1][len(costArray[0]) - 1])


		

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = costArray[len(costArray) - 1][len(costArray[0]) - 1];
		alignment1 = alignSeq1
		alignment2 = alignSeq2
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
	

	def align_banded(self, seq1, seq2, align_length):
		if (abs(len(seq1) - len(seq2)) > align_length):
			#can't find the cost for this band length
			score = math.inf;
			alignment1 = 'No Alignment Possible'
			alignment2 = 'No Alignment Possible'
			return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}



		#print('starting align ************************************')
		#print(seq1)
		#print(seq2)
		left = 1
		up = 2
		diagonal = 3
		delin = 5
		sub = 1
		found = -3

		costArray = []
		prevArray = []
		array = []
		for i in range(align_length + align_length - 1):
			array.append(0)
		for j in range(len(seq2) + 1):
			costArray.append(array.copy())
			prevArray.append(array.copy())

		for i in range(align_length):
			costArray[0][i + align_length - 1] = i*delin
			if (i != 0):
				prevArray[0][i + align_length - 1] = left
		for j in range(min(align_length, len(seq2) + 1)):
			costArray[j][align_length - j - 1] = j*delin
			if (j != 0):
				prevArray[j][align_length - j - 1] = up
		#print(costArray)
		#print('filling out tables----------------------------------')
		for j in range(1, len(seq2) + 1):
			beginInd = max((j - align_length + 1), 1)
			endInd = min(beginInd + align_length + j - 1, beginInd + align_length + align_length - 1, len(seq1) + 1)
			for i in range(beginInd, endInd):
				#debugging
				#print(str(j) + ' ' + str(i) + ' (' + seq2[j - 1] + ' ' + seq1[i - 1] + ')')
				costColInd = max(align_length - j + i - 1, 0)
				#print(costColInd)
				subCost = costArray[j - 1][costColInd] + sub
				if (seq1[i - 1] == seq2[j - 1]):
					subCost = subCost - sub + found
				
				#assign in opposite order of wanted left, up, diagonal
				minimum = subCost
				prev = diagonal
				if (costColInd + 1 < len(costArray[0])):
					cost = costArray[j-1][costColInd + 1] + delin
					if (minimum > cost):
						minimum = cost
						prev = up
				if (costColInd - 1 >= 0):
					cost = costArray[j][costColInd - 1] + delin
					if (minimum > cost):
						minimum = cost
						prev = left
				
				costArray[j][costColInd] = minimum
				prevArray[j][costColInd] = prev
		
		#print('filled out arrays-----------------')
		#for j in range(len(costArray)):
		#	print(costArray[j])
		#print('prev')
		#for j in range(len(prevArray)):
		#	print(prevArray[j])
		
		cost = math.inf
		i = -1
		j = len(seq2)
		if (len(seq2) - len(seq1) >= 0):
			i = align_length - (len(seq2) - len(seq1))
			cost = costArray[j][i]
		else:
			i = align_length
			cost = costArray[j][i]
		#print(cost)

		#getting the alignment
		alignSeq1 = ''
		alignSeq2 = ''
		#print(str(j) + ' ' + str(i))
		#print(str(j - 1) + ' ' + str(i + j - align_length))
		#print()
		while (prevArray[j][i] != 0):
			if (prevArray[j][i] == left):
				alignSeq1 = seq1[i + j - align_length] + alignSeq1
				alignSeq2 = '-' + alignSeq2
				i -= 1
			elif (prevArray[j][i] == up):
				alignSeq1 = '-' + alignSeq1
				alignSeq2 = seq2[j-1] + alignSeq2
				i += 1
				j -= 1
			elif (prevArray[j][i] == diagonal):
				alignSeq1 = seq1[i + j - align_length] + alignSeq1
				alignSeq2 = seq2[j-1] + alignSeq2
				j -= 1

			#debugging
			#print(str(j) + ' ' + str(i))
			#print(str(j - 1) + ' ' + str(i + j - align_length))
			#print()

			#trim so that only the first 100 characters are visible
			if (len(alignSeq1) > 100):
				alignSeq1 = alignSeq1[0:len(alignSeq1) - 1]
			if (len(alignSeq2) > 100):
				alignSeq2 = alignSeq2[0:len(alignSeq2) - 1]

		print('results*********')
		print(cost)
		print(alignSeq1)
		print(alignSeq2)

		score = cost
		alignment1 = alignSeq1
		alignment2 = alignSeq2
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}


