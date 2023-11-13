#!/usr/bin/python
# coding=utf-8
#distutils: language = c++
#cython: language_level=3

import cython
import numpy as np 
cimport cython 
cimport numpy as np
from cython.parallel cimport prange
from cython.parallel cimport parallel


from libc.string cimport memset
from libc.stdlib cimport malloc as Malloc, free as Free, realloc as Realloc
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free



DTYPE=np.int32
ctypedef np.int32_t DTYPE_t
ctypedef unsigned char uchar
ctypedef unsigned short int16

blocksize = 20


cdef class LinkArray:
	
	cdef int16** data
	cdef int size
	
	def __cinit__(self):
		
		self.size = 0
		self.data = <int16**> Malloc(sizeof(int16*))
		#self.data[0] = <int16*> Malloc(blocksize*sizeof(int16))
		
		
	def __dealloc__(self):
		
		for blockindex in xrange(0, self.size//blocksize):
			
			Free(self.data[blockindex])
			
	def add(self, int16 element):
		
		blocknum = self.size//blocksize
		
		if blocknum % blocksize == 0 :
			
			self.data = <int16**> Realloc(self.data,(blocknum+blocksize)*sizeof(int16*))
			
		if  self.size % blocksize == 0:
			
			self.data[blocknum] = <int16*> Malloc(blocksize*sizeof(int16))
			
		self.data[blocknum][self.size%blocksize] = element
		
		self.size += 1
		
	def get(self, int index):
		
		return self.data[index//blocksize][index % blocksize]

	def getsize(self):
		return self.size
	
	def getrow(self, np.ndarray[int, ndim=1, mode = 'c'] keys):
		cdef int16* block
		cdef int16 keyssize = 0 
		cdef int16 key
                
		for blockindex in xrange(0, self.size//blocksize):
                        
			block = self.data[blockindex]
                        
			for inblockindex in xrange(blocksize):
                                
				keys[keyssize] = block[inblockindex]
				keyssize += 1                   
             
		block = self.data[self.size//blocksize]
                
		for inblockindex in xrange(self.size%blocksize):
                        
			keys[keyssize] = block[inblockindex]
			keyssize += 1 
                                
		return keyssize				

	
	def count(self, np.ndarray[int, ndim=1, mode = 'c'] keys, np.ndarray[int, ndim=1, mode = 'c'] counter, set typeset):
		
		cdef int16* block
		cdef int16 keyssize = 0 
		cdef int16 key
		
		for blockindex in xrange(0, self.size//blocksize):
			
			block = self.data[blockindex]
			
			for inblockindex in xrange(blocksize):
				
				key = block[inblockindex]
				
				if key in typeset:
					
					counter[key] += 1
					
				else:
					
					typeset.add(key)
					counter[key] = 1
					
					keys[keyssize] = key
					keyssize += 1
					
		block = self.data[self.size//blocksize]
		
		for inblockindex in xrange(self.size%blocksize):
			
			key = block[inblockindex]
			
			if key in typeset:
				
				counter[key] += 1
				
			else:
				
				typeset.add(key)
				counter[key] = 1
				
				keys[keyssize] = key
				keyssize += 1
				
				
		return keyssize
	
cdef void findmatch_L2(int gnamesize, float querycounts, int* keys, int keysize, int* counter, float* match_counts):
	
	cdef int count_i, index_start, key_i, key_j, i , j
	
	cdef float delta
	
	delta = querycounts

	for i in range(0,keysize):
		
		key_i = keys[i]
	
		count_i = counter[key_i]

		index_start = gnamesize*key_i
		
		match_counts[index_start + key_i] += 0.5*float(count_i * count_i)/delta
		
		for j in range(i+1,keysize):
			
			key_j = keys[j]
			
			match_counts[index_start + key_j] += float(count_i * counter[key_j])/delta


cdef void findmatch_L1(int gnamesize, float querycounts, int* keys, int keysize, int* counter, float* match_counts) nogil:
	
	cdef int count_i, index_start, key_i, key_j, i , j
	
	cdef float delta
	
	delta = querycounts
	
	for i in range(0,keysize):
		
		key_i = keys[i]
		
		count_i = counter[key_i]
		
		index_start = gnamesize*key_i
		
		match_counts[index_start + key_i] += count_i
		
		for j in range(i+1,keysize):
			
			key_j = keys[j]
			
			match_counts[index_start + key_j] += min(count_i, counter[key_j])
			
cdef void findmatch_L1_LM(int gnamesize, float querycounts, int* keys, int keysize, int* counter, float* match_counts) nogil:
	
	cdef int count_i, index_start, key_i, key_j, i , j
	
	cdef float delta
	
	delta = querycounts
	
	for i in range(0,keysize):
		
		key_i = keys[i]
		
		count_i = counter[key_i]
		
		index_start = gnamesize*key_i - (key_i*(key_i+1)//2)
		
		match_counts[index_start + key_i] += count_i
		
		for j in range(i+1,keysize):
			
			key_j = keys[j]
			
			match_counts[index_start + key_j] += min(count_i, counter[key_j])


			
cdef class SparseKmerMartrix:
	
	cdef int kmersize, gnamesize, querysize
	kmerrows = []
	querycounts = {}
	
	def __cinit__(self, int gnamesize , int kmersize = 1000, int querysize = -1):
		
		self.kmersize = kmersize 
		self.querysize = gnamesize if querysize == -1 else querysize
		self.gnamesize = gnamesize
		
		
		for index in xrange(0,self.kmersize):
			
			self.kmerrows.append(LinkArray())
			
	def __dealloc__(self):
		
		for index in xrange(0,self.kmersize):
			
			del self.kmerrows[0]
			
	def addkmer(self):
		
		self.kmerrows.append(LinkArray())
		
		self.kmersize += 1
		
	def addquerycounts(self, kindex, count = 1.0):
		
		self.querycounts[kindex] = self.querycounts.get(kindex, 0) + count

	def getkmerrow(self, int kindex):

		cdef np.ndarray[int, ndim=1, mode = 'c'] counter = np.zeros(self.kmerrows[kindex].getsize(), dtype=np.int32)
                
		self.kmerrows[kindex].getrow(counter)

		return counter
			
	
	def getgenenumber(self):
		
		return self.gnamesize
	
	def add(self, int kmerindex, int16 gnameindex):
		
		self.kmerrows[kmerindex].add(gnameindex)
		
	def getkmercounts(self, int kindex):
		
		cdef np.ndarray[int, ndim=1, mode = 'c'] counter = np.zeros(self.gnamesize, dtype=np.int32)
		
		cdef np.ndarray[int, ndim=1, mode = 'c'] keys = np.zeros(self.gnamesize, dtype=np.int32)
		
		typeset = set({})
		
		keysize = self.kmerrows[kindex].count(keys, counter,typeset)
		
		return counter

	def getkmersamples(self, int kindex):
	
		cdef int index
		cdef int sampleindex
		cdef int size 
		
		size = self.kmerrows[kindex].getsize()	
	
		samplecounts = {}	
		for index in xrange(size):
			
			sampleindex = self.kmerrows[kindex].get(index)
			
			if sampleindex in samplecounts:
				
				samplecounts[sampleindex] += 1
				
			else:
				
				samplecounts[sampleindex] = 1
				
		return samplecounts	
	
	def SquareMatrix(self, weighted = 0):
		
		cdef np.ndarray[float, ndim=1, mode = 'c'] match_counts =  np.zeros(self.gnamesize * self.gnamesize, dtype=np.float32)
		cdef float* match_counts_prt = <float*> match_counts.data
		
		cdef np.ndarray[int, ndim=1, mode = 'c'] counter = np.zeros(self.gnamesize, dtype=np.int32)
		cdef int* counter_prt = <int*> counter.data
		
		cdef np.ndarray[int, ndim=1, mode = 'c'] keys = np.zeros(self.gnamesize, dtype=np.int32)
		cdef int* keys_prt = <int*> keys.data
		cdef int keysize = 0
		cdef int kmersize = self.kmersize
		cdef int kindex
		cdef float querycount,default
		cdef int count_i, index_start, key_i, key_j, i , j
		
		typeset = set({})
	
		for kindex in prange(kmersize, schedule='dynamic', nogil=1):
			
			with gil:
				
				typeset.clear()
				
				keysize = self.kmerrows[kindex].count(keys, counter,typeset)
				
				default = 1.0
				if len(typeset) < 2 :
					default = 20.0
					
				if weighted:
					querycount = self.querycounts.get(kindex, default)
				else:
					querycount = default
					
				if querycount:
					findmatch_L2(self.gnamesize, querycount, keys_prt, keysize, counter_prt, match_counts_prt) 
					
		for i in range(0,self.gnamesize):
			
			index_start = self.gnamesize*i
			
			match_counts[index_start + i] = 2*match_counts[index_start + i] 
			
			for j in range(i+1,self.gnamesize):
				
				match_counts[j*self.gnamesize + i] = match_counts[index_start + j]
				
				
		return match_counts
	
	
	def MatchMatrix(self):
		
		cdef np.ndarray[float, ndim=1, mode = 'c'] match_counts =  np.zeros(self.gnamesize * self.gnamesize, dtype=np.float32)
		cdef float* match_counts_prt = <float*> match_counts.data
		
		cdef np.ndarray[int, ndim=1, mode = 'c'] counter = np.zeros(self.gnamesize, dtype=np.int32)
		cdef int* counter_prt = <int*> counter.data
		
		cdef np.ndarray[int, ndim=1, mode = 'c'] keys = np.zeros(self.gnamesize, dtype=np.int32)
		cdef int* keys_prt = <int*> keys.data
		cdef int keysize = 0
		cdef int kmersize = self.kmersize
		cdef int kindex
		cdef float querycount
		cdef int count_i, index_start, key_i, key_j, i , j
		
		typeset = set({})
		
		for kindex in prange(kmersize, schedule='dynamic', nogil=1):
			
			with gil:
				
				typeset.clear()
				
				keysize = self.kmerrows[kindex].count(keys, counter,typeset)
				
				if len(typeset) < 2:
					continue			
	
				findmatch_L1(self.gnamesize, self.querycounts.get(kindex, 1.0), keys_prt, keysize, counter_prt, match_counts_prt) 
				
				
		for i in range(0,self.gnamesize):
			
			index_start = self.gnamesize*i
						
			for j in range(i+1,self.gnamesize):
				
				match_counts[j*self.gnamesize + i] = match_counts[index_start + j]
				
				
		return match_counts

	def MatchMatrix_LM(self):
		
		cdef np.ndarray[float, ndim=1, mode = 'c'] match_counts =  np.zeros(self.gnamesize * self.gnamesize - (self.gnamesize * (self.gnamesize -1)) // 2 , dtype=np.float32)
		cdef float* match_counts_prt = <float*> match_counts.data
		
		cdef np.ndarray[int, ndim=1, mode = 'c'] counter = np.zeros(self.gnamesize, dtype=np.int32)
		cdef int* counter_prt = <int*> counter.data
		
		cdef np.ndarray[int, ndim=1, mode = 'c'] keys = np.zeros(self.gnamesize, dtype=np.int32)
		cdef int* keys_prt = <int*> keys.data
		cdef int keysize = 0
		cdef int kmersize = self.kmersize
		cdef int kindex
		cdef float querycount
		cdef int count_i, index_start, key_i, key_j, i , j
		
		typeset = set({})
		
		for kindex in prange(kmersize, schedule='dynamic', nogil=1):
			
			with gil:
				
				typeset.clear()
				
				keysize = self.kmerrows[kindex].count(keys, counter,typeset)
				
				if len(typeset) < 2:
					continue			
				
				findmatch_L1_LM(self.gnamesize, self.querycounts.get(kindex, 1.0), keys_prt, keysize, counter_prt, match_counts_prt) 
				
				
				
		return match_counts	
