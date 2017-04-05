import numpy
n_qtl_interval_genes = 210
n_multihit_genes = 27

n_yeast_genes = 5858

n_iter = 100000
n_overlap = []

for i in range(n_iter):
	
	choice1 = numpy.random.choice(numpy.arange(n_yeast_genes), n_multihit_genes)
	choice2 = numpy.random.choice(numpy.arange(n_yeast_genes), n_qtl_interval_genes)
	
	n_overlap.append( len( set(choice1) & set(choice2) ) )
n_overlap = numpy.array(n_overlap)	
print numpy.sum(n_overlap > .5)/float(n_iter)
print numpy.sum(n_overlap > 1.5)/float(n_iter)
print numpy.sum(n_overlap > 2.5)/float(n_iter)

print numpy.sum(n_overlap)/float(n_iter)