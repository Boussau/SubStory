Substory needs to be run after bppml. An example is provided in the Example folder.

##########################################
## Meaning of the options to the program:
##########################################

## Ancestral sequence reconstruction method (no other option for the moment):
asr.method=marginal 

## Report probabilities:
asr.probabilities=yes 

## If we want to sample from the posterior distribution of ancestral sequences, and get N=10 sequences per node, use these two options: 
asr.sample=yes
# How many samples should we use?
asr.sample.number=10

## Shall we add extant sequences to output file?
output.nodes.add_extant = yes

## Alignment information log file (site specific rates, probabilities, etc):
output.sites.file = $(DATA).sites.csv

## Nodes information log file
output.nodes.file = $(DATA).nodes.csv

## Write sequences:
output.sequence.file = $(DATA).ancestors.fasta
output.sequence.format = Fasta

## What type of substitutions do we want to keep track of?
map.type=All

## Output all substitutions per branch (lines) and per type (columns)
output.mapping.file=substitutions.csv
