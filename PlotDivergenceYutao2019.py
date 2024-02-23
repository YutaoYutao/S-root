import os

dirname = '/Users/zyt/Desktop/'
fastaname = 'AccessionLocus20191220.fasta'

def replace_dots_with_Xes(infilename, outfilename):
	
	id2seq, idlist = fasta2seqhash_and_idlist(open(infilename))
	outfile = open(outfilename, 'w')
	for id in idlist:
		outfile.write('>'+id+'\n'+id2seq[id].upper().replace('.', 'X')+'\n')
	outfile.close()
	
# read fasta and store it in a dictionary linking headers to sequence
# keep the order of sequences in a list (order is not kept in a dictionary)
def fasta2seqhash_and_idlist(fastafile):
	
	line        = fastafile.readline()
	currentgene = currentgene = line[1:]
	
	seq=''
	id2seq = {}
	line     = fastafile.readline()
	idlist = []
	while len(line)>0:
		if line[0]=='>':
			if id2seq.has_key(currentgene) and id2seq[currentgene] != seq: print 'key', currentgene, "already exists, with a different sequence! I will overwrite it. I suggest you will set keep_whole_header to 'True'"
			
			id2seq[currentgene]=seq
			idlist.append(currentgene)
			#intialize new
			currentgene = line[1:]
			seq=''
		else:
			#add this line to sequence
			seq += line.strip()
		line = fastafile.readline()
	
	#don't forget the last one :-)
	id2seq[currentgene]=seq
	idlist.append(currentgene)
	
	return id2seq, idlist
	
	
	

def divergence_and_gaps_per_position(id2seq, idlist):
	
	refid  = idlist[0]
	refseq = id2seq[refid].upper()
	
	#for each position, keep a list where you keep divergence of all other sequences
	similarity = []
	Ngaps      = []
	Nmissing   = []
	totalpos   = []
	
	for pos in refseq:
		if pos != '-':
			similarity.append(0)
			Ngaps.append(0)
			Nmissing.append(0)
			totalpos.append(0)
	
	# now go through all sequences	
	#print refid, refseq
	for id in idlist[1:]: #[1:] means skip the first one in the list (because this is the reference)
		
		seq = id2seq[id].upper()
		#print id, seq
		
		seqindex = 0
		refindex = 0
		for refnuc in refseq:
			if refnuc != '-':
				#print seq[i], refnuc
				if seq[seqindex] != 'X': #if a position is missing, just ignore the whole thing in divergence stats
					if seq[seqindex] == refnuc: 
						similarity[refindex]+=1
						
					elif seq[seqindex] == '-' : 
						Ngaps[refindex] += 1
					totalpos[refindex] += 1
					
				else: Nmissing[refindex] += 1 #count how often data is missing for a certain position
				refindex += 1
				
			seqindex+=1
			
	
	prop_same_per_position    = []
	prop_gaps_per_position    = []
	prop_missing_per_position = []
	for i in range(len(refseq.replace('-', ''))):
		prop_same_per_position.append(similarity[i]/float(totalpos[i]))
		prop_gaps_per_position.append(Ngaps[i]/float(totalpos[i]))
		prop_missing_per_position.append(Nmissing[i]/float(totalpos[i]))
	
	return prop_same_per_position, prop_gaps_per_position, prop_missing_per_position
	
	
	

def test():
	
	id2seq = {}
	
	id2seq['r'] = 'abcdEFG----hijkmlnopQSRTUVWxyz'
	id2seq['1'] = 'abcdEFG----hijkmlnopQSRTUVWxyz'
	id2seq['2'] = 'abcdEFG----hijkmlnopQSRTUVWxyz'
	id2seq['3'] = 'abcdEFG----hijkmlnopQSRTUVWxyz'
	idlist = ['r', '1', '2', '3']
	
	print "All the same, none missing:"
	print divergence_and_gaps_per_position(id2seq, idlist)
	
	id2seq = {}
	id2seq['r'] = 'abcdEFG----hijkmlnopQSRTUVWxyz'
	id2seq['1'] = 'Ab--EFG----hijkmlnopQSRTUVWxyz'
	id2seq['2'] = 'ab--EFGgggghijkmlnopQSRTUVWxyz'
	id2seq['3'] = 'abcdEFGgghhhijkmlnopQSR--VWxyz'
	
	print "Some difference, none missing:"
	print divergence_and_gaps_per_position(id2seq, idlist)
	
	id2seq = {}
	id2seq['r'] = 'abcdEFG----hijkmlnopQSRTUVWxyz'
	id2seq['1'] = 'Ab--EFG----hijkmlnXXQSRTUVWxyz'
	id2seq['2'] = 'ab--EFGgggghijkmlnopQSRTUVWxyz'
	id2seq['3'] = 'XXXXXFGgghhhijkmlnopQSR--VWxyz'
	
	print "Some difference, some missing:"
	print divergence_and_gaps_per_position(id2seq, idlist)
	
	
	
	
def plot(outfile_base, prop_same_per_position, prop_gaps_per_position, prop_missing_per_position, orf2exonpositions):
	
	gnufilename  = outfile_base+'.multiplot.gnu'
	
	# write data on similarites etc per position to a file, so gnuplot can read it
	datafilename = outfile_base+'.dat'
	datafile     = open(datafilename, 'w')
	for i in range(len(prop_same_per_position)):
		datafile.write(str(i+1)+'\t'+str(prop_same_per_position[i])+'\t'+str(prop_gaps_per_position[i])+'\t'+str(prop_missing_per_position[i])+'\n')
	datafile.close()
	
	# write data on similarites etc per position to a file, so gnuplot can read it
	labelstr = ''
	
		
	for orf in orf2exonpositions.keys():
		
		orffilename = outfile_base+orf+'.dat'
		orffile     = open(orffilename, 'w')
		for (start,end) in orf2exonpositions[orf]:
			orffile.write(str(start)+'\t1\n'+str(end)+'\t1\n\n')
			
		labelstr += 'set label "'+orf+'" at '+str(orf2exonpositions[orf][0][0]+60) + ', 1.3\n'
		
		orffile.close()
	
	gnu = "set terminal svg\nset output '"+outfile_base+".svg'\n\n"
	gnu += '# uncomment the lines below and comment out the lines above to get an eps plot\n'
	gnu += "# set terminal postscript eps enhanced color\n# set output '"+outfile_base+".eps'\n\n"
	#
	# Set top and bottom margins to 0 so that there is no space between plots.
	# Fix left and right margins to make sure that the alignment is perfect.
	# Turn off xtics for all plots except the bottom one.
	# In order to leave room for axis and tic labels underneath, we ask for
	# a 5-plot layout but only use the top 4 slots.
	#
	gnu += 'set tmargin 1\n'
	gnu += 'set bmargin 0\n'
	gnu += 'set lmargin 8\n'
	gnu += 'set rmargin 3\n'
	gnu += 'set tics scale 0 font ",8"\n'
	gnu += 'unset key\n'
	gnu += 'unset xtics\n'
	gnu += 'set multiplot layout 5,1\n'
	
	gnu += 'set yrange [-0.1:1.1]\n'
	gnu += 'set ylabel "missing data"\n'
	
	gnu += "plot '"+datafilename+"' using 1:4 w l lt 1\n" # missing data
	gnu += 'set ylabel "gaps"\n'
	gnu += "plot '"+datafilename+"' using 1:3 w l lt 2\n" # gaps
	gnu += 'set ylabel "similarity"\n'
	gnu += "plot '"+datafilename+"' using 1:2 w l lt 3\n" # similarity
	
	nlinetypes = 7
	gnu += 'set yrange [0:2]\n'
	gnu += 'set style line 1 lt 1 lw 5 pt 7 ps 0.2 lc rgb "blue"\n'
	gnu += 'set style line 2 lt 1 lw 5 pt 7 ps 0.2 lc rgb "dark-green"\n'
	gnu += 'set style line 3 lt 1 lw 5 pt 7 ps 0.2 lc rgb "orange"\n'
	gnu += 'set style line 4 lt 1 lw 5 pt 7 ps 0.2 lc rgb "purple"\n'
	gnu += 'set style line 5 lt 1 lw 5 pt 7 ps 0.2 lc rgb "brown"\n'
	gnu += 'set style line 6 lt 1 lw 5 pt 7 ps 0.2 lc rgb "black"\n'
	gnu += 'set style line 7 lt 1 lw 5 pt 7 ps 0.2 lc rgb "red"\n'

	gnu += 'set xtics nomirror\nunset ytics\n'
	gnu += labelstr
	gnu += "plot '"
	for orfindex, orf in enumerate(orf2exonpositions.keys()):
		orffilename = outfile_base+orf+'.dat'
		ls = orfindex + 1
		if ls > nlinetypes: ls = ls%nlinetypes
		
		gnu += orffilename+"' w l ls "+str(ls)+", '"
		
	gnu = gnu[:-3]
		
	gnufile = open(gnufilename, 'w')
	gnufile.write(gnu)
	gnufile.close()
	return gnufilename
	




def make_plot1():
	fastafilename  = dirname + 'Locus2.aligned.clustalo.fasta'
	id2seq, idlist = fasta2seqhash_and_idlist(open(fastafilename))
	prop_same_per_position, prop_gaps_per_position, prop_missing_per_position = divergence_and_gaps_per_position(id2seq, idlist)
	
	outfile_base = dirname + 'LOCUS2_PLOT'
	orf2exonpositions = {}
	orf2exonpositions['AT3g57620'] = [(565,2208)]
	orf2exonpositions['AT3g57630'] = [(2544,2834),(2918,3016),(3125,3356),(3468,3548),(3644,3891),(3990,4103),(4188,4322),(4428,4514),(4613,4828),(4937,5192),(5283,5469),(5568,5661),(5740,6081)]
	orf2exonpositions['AT3g57640'] = [(7922,8992)]
	orf2exonpositions['AT3g57645'] = [(9343,10205)]
	orf2exonpositions['AT3g57650'] = [(12752,12829),(13197,13319),(13429,13468),(13583,13656),(13952,14011),(14232,14309),(14470,14628),(14835,14960),(15204,15296),(15429,15515),(15589,15840)]
	
		
	plot(outfile_base, prop_same_per_position, prop_gaps_per_position, prop_missing_per_position, orf2exonpositions)
	
	
def make_plot2():
	infilename  = dirname + fastaname
	outfilename = dirname + fastaname.replace('.fasta', '.dots2x.fasta')
	replace_dots_with_Xes(infilename, outfilename)
	
	aligned_filename = outfilename.replace('.fasta', '.aligned.clustalo.fasta')
	cmnd = '/usr/local/bin/clustalo -i '+outfilename+' -t DNA > '+aligned_filename
	#cmnd = 'sudo chgrp -R admin /usr/local/bin/clustalo -i '+outfilename+' -t DNA > '+aligned_filename
	print cmnd
	print os.system(cmnd)
	
	id2seq, idlist = fasta2seqhash_and_idlist(open(aligned_filename))
	prop_same_per_position, prop_gaps_per_position, prop_missing_per_position = divergence_and_gaps_per_position(id2seq, idlist)
	
	outfile_base = outfilename.split('.aligned.clustalo.fasta')[0]
	
	orf2exonpositions = {}
	orf2exonpositions['AT3g57620'] = [(565,2208)]
	orf2exonpositions['AT3g57630'] = [(2544,2834),(2918,3016),(3125,3356),(3468,3548),(3644,3891),(3990,4103),(4188,4322),(4428,4514),(4613,4828),(4937,5192),(5283,5469),(5568,5661),(5740,6081)]
	orf2exonpositions['AT3g57640'] = [(7922,8992)]
	orf2exonpositions['AT3g57645'] = [(9343,10205)]
	orf2exonpositions['AT3g57650'] = [(12752,12829),(13197,13319),(13429,13468),(13583,13656),(13952,14011),(14232,14309),(14470,14628),(14835,14960),(15204,15296),(15429,15515),(15589,15840)]
		
	gnufilename = plot(outfile_base, prop_same_per_position, prop_gaps_per_position, prop_missing_per_position, orf2exonpositions)
	
	cmnd = 'gnuplot '+gnufilename
	print cmnd, os.system(cmnd)
	
	
if __name__ == "__main__":
	#test()
	#make_plot1()
	
	make_plot2()
