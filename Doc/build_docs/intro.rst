.. intro:

**********************************
Motivation for profiling DGE tools
**********************************

"Which of my genes are differentially expressed?" is one of the most commonly 
asked questions in modern cell biology. For RNA-seq the DGE work-flow typically
involves extracting and sequencing the RNA from the experimental 
conditions, normalizing and cleaning the data, calculating a fold-change for 
each gene between the conditions in question, and assigning this differential 
expression a statistical significance based on some assumed properties of the 
data. This work-flow is, of course, subject to errors and biases at all stages 
and scientists put a lot of effort into trying to understand, control and 
minimise these effects. 

The final step of this work-flow, calculating the fold-change for each gene, 
and assigning a statistical significance to the 
fold-change, has often been assumed to be a solved problem and has not yet 
received the same measure of scrutiny and as other steps in this process. 
Numerous methods exist for performing this final step but unfortunately there 
is no strong consensus in the community about what is the most appropriate 
method to use to calculate DGE from RNA-seq data or, indeed, whether any of 
them produce consistent, reproducible, accurate results at all. Worse, typical 
papers that use RNA-Seq data usually only present differential expression 
results from one of these methods, without any estimation of how these results
depend on the method used to calculate the differential expression.

This experiment was designed and performed specifically to address some of these
issues.

..	toctree::
	:maxdepth: 1

	install