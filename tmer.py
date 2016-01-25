import pickle
import scipy
import numpy as np
import sys
from scipy.sparse import dok_matrix
KMER_LEN=50

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

def erase():
	print(CURSOR_UP_ONE + ERASE_LINE)

def enumerateKmers(text):
	kmers = []
	for i in range(0,len(text) - KMER_LEN):
		kmers.append(text[i:(i+KMER_LEN)])
	return kmers

class TranscriptClassifier:
	tr_ids = []
	tr_seqs = []
	kmer_catalog = {}
	source_filename=None
	adj=None

	def loadFasta(self,filename):
		self.source_filename = filename
		lines = open(filename,"r").read().split("\n")

		curr_seq=False
		for line in lines:
			l = line.strip()
			if len(l) == 0:
				continue

			if l[0]==">":
				self.tr_ids.append(l)
				if curr_seq != False:
					self.tr_seqs.append(curr_seq)
				curr_seq = ""
			else:
				curr_seq += l

		if curr_seq != "":
			self.tr_seqs.append(curr_seq)

		print("Loaded %d sequences from %s"%(len(self.tr_seqs),filename))

	def buildKmerCatalog(self):
		for i in range(len(self.tr_seqs)):
			print("Build catalog for %d"%i)
			#erase()
			for kmer in enumerateKmers(self.tr_seqs[i]):
				if not kmer in self.kmer_catalog:
					self.kmer_catalog[kmer]=set()
				self.kmer_catalog[kmer].add(i)

		print("Built kmer catalog (%d entries)"%len(self.kmer_catalog))
		#pickle.dump(self.kmer_catalog,open(self.source_filename+"-kmers.dat","w"))


	def getMates(self,i):
		kmers = enumerateKmers(self.tr_seqs[i])
		mates = set()

		for kmer in kmers:
			#print("Searching kmer %s"%kmer)
			hits = self.kmer_catalog[kmer]
			mates = set(list(mates)+list(hits))

		return list(mates)

	def filterByLength(self,min_length):
		new_tr_ids=[]
		new_tr_seqs=[]
		for i in range(len(self.tr_seqs)):
			if len(self.tr_seqs[i]) >= min_length:
				new_tr_ids.append(self.tr_ids[i])
				new_tr_seqs.append(self.tr_seqs[i])
		
		self.tr_ids=new_tr_ids
		self.tr_seqs=new_tr_seqs

	def write(self,filename):
		handle = open(filename,"w")
		for i in range(len(self.tr_seqs)):
			handle.write(self.tr_ids[i])

			for j in range(len(self.tr_seqs[i])):
				if j % 60 == 0:
					handle.write("\n")
				handle.write(self.tr_seqs[i][j])
			handle.write("\n")
		handle.close()

	def buildAdjMatrix(self):
		print("Building adj. matrix")
		#self.adj = [[0 for j in range(len(self.tr_seqs))] for i in range(len(self.tr_seqs))]
		self.adj = dok_matrix((len(self.tr_seqs),len(self.tr_seqs)),dtype=np.float32)
		self.adj_steps=[]
		self.not_null_prev = 0
		count=0
		for i in range(len(self.tr_seqs)):
			mates=self.getMates(i)
			for m in mates:
				self.adj[i,m]=1
				count+=1

		self.initial_connection = self.adj.getnnz()

		self.adj_steps.append(self.adj.copy())
		self.initmatrix = self.adj.copy()
		self.final_matrix = self.adj.copy()
		print("Finished building adj. matrix, connections: %d (check %d)"%(self.initial_connection,count))

	def stepAdjMatrix(self):
		print("Computing adj. matrix metric")
		not_null=self.final_matrix.getnnz()
		#for i in range(len(self.tr_seqs)):
		#	for j in range(len(self.tr_seqs)):
		#		if self.adj[i,j]!=0:
		#			not_null+=1

		print("Elements in adj. matrix: %d, previous: %d"%(not_null,self.not_null_prev))
		if not_null == self.not_null_prev:
			#print("Computing reachability matrix from steps")
			self.adj = self.final_matrix
			#self.adj = dok_matrix((len(self.tr_seqs),len(self.tr_seqs)),dtype=np.float32)
			#for m in self.adj_steps:
			#	self.adj = self.adj + m
			self.final_connection = self.adj.getnnz()
			print("Initial connection: %d, Final connection: %d"%(self.initial_connection,self.final_connection))
			
			#print("Checking matrix")
			#for i in range(len(self.tr_seqs)):
			#	for j in range(len(self.tr_seqs)):
			#		if self.initmatrix[i,j] > 0 and (not self.adj[i,j]>0):
			#			print("Init > 0 but new not > 0!")
			#			throw

			return False

		self.not_null_prev = not_null

		print("Computing power of adj. matrix")
		self.adj=self.adj**2
		self.adj_steps.append(self.adj.copy())
		self.final_matrix = self.final_matrix + self.adj
		print("Finished computing power of adj. matrix")
		print("")
		return True

	def getMatesAdjMatrix(self,i):
		mates=[]
		for j in range(len(self.tr_seqs)):
			if self.adj[i,j] > 0:
				mates.append(j)
		return mates

	def getGroupIdentifier(self,i):
		rows,cols = self.adj.getrow(i).nonzero()
		gid = ""
		gc=0
		j=0
		for col in cols:
			gid+=str(col)+"|"
			if col == i:
				gc=j
			j+=1
		return gid,gc,j	

	def writeAnnotated(self,filename):
		print("Re-annotating sequences with group identifiers (fastmode)...")
		handle = open(filename,"w")
		handle_conv = open(filename+".conv","w")
		handle_gsizes=open(filename+".groupsizes","w")
		handle_gsizedist=open(filename+".groupsizes.dist","w")

		group_sizes={}
		groups={}
		gts={}
		group_identifier = 0
		for i in range(len(self.tr_seqs)):
			gid,gc,gt=self.getGroupIdentifier(i)
			
			if gt > 1:
				if gid in groups:
					ac_gid = groups[gid]
				else:
					groups[gid]=group_identifier
					ac_gid = group_identifier
					group_identifier += 1
					if not gt in gts:
						gts[gt]=0
					gts[gt]+=1
			else:
				ac_gid = group_identifier
				group_identifier += 1

			ac_gid += 1
			ac_sgid = gc+1
			text_id="C%d_%d"%(ac_gid,ac_sgid)
			print(text_id)
			handle.write(">%s"%text_id)
			handle_conv.write("%s,%s\n"%(self.tr_ids[i].split(" ")[0][1:],text_id))
			handle_gsizes.write("C%d,%d\n"%(ac_gid,gt))

                        for j in range(len(self.tr_seqs[i])):
                                if j % 60 == 0:
                                        handle.write("\n")
                                handle.write(self.tr_seqs[i][j])
                        handle.write("\n")

		handle_gsizedist.write("Done. Number of clusters: %d\n"%len(groups))
		handle_gsizedist.write("Cluster Size | # Clusters\n")
		for k in gts:
			handle_gsizedist.write("%4d | %4d\n"%(k,gts[k]))

		handle.close()
		handle_conv.close()
		handle_gsizes.close()		

	def writeAnnotatedOld(self,filename):
		print("Re-annotating sequences with group identifiers...")
		group_id = 1
		groups = {}
		group_counters = {}

		handle = open(filename,"w")
		handle_conv = open(filename+".conv","w")
		for i in range(len(self.tr_seqs)):
			mates = self.getMatesAdjMatrix(i)

			if i in groups:
				curr_id = groups[i]
				curr_sub_id = group_counters[curr_id]
				group_counters[curr_id] = group_counters[curr_id] + 1
			else:
				for m in mates:
					groups[m]=group_id
				groups[i]=group_id
				
				if group_id in group_counters:
					throw

				group_counters[group_id]=2
				curr_id=group_id
				curr_sub_id=1

				group_id += 1

			#handle.write("%s%s"%(self.tr_ids[i],str(mates)))
			text_id = "G%d_%d"%(curr_id,curr_sub_id)
			handle.write(">%s"%text_id)
			handle_conv.write("%s,%s\n"%(self.tr_ids[i].split(" ")[0][1:],text_id))

                        for j in range(len(self.tr_seqs[i])):
                                if j % 60 == 0:
                                        handle.write("\n")
                                handle.write(self.tr_seqs[i][j])
                        handle.write("\n")

			#print(text_id)
			#erase()
		handle.close()
		handle_conv.close()

if len(sys.argv) < 3:
	print("Usage: tmer.py <Input Fasta> <Output Clustered Fasta>")

main = TranscriptClassifier()
main.loadFasta(sys.argv[1])
main.buildKmerCatalog()

main.buildAdjMatrix()


while main.stepAdjMatrix():
	pass

main.writeAnnotated(sys.argv[2])

