# ASRUs_generator

ASRU stands for Alternatively Spliced Repetitive Unit. An ASRU is a connected component in the similarity subgraph of an Evolutionnary Splicing Graph (=ESG, [1]). 
For a given gene we associate to it two graphs, an ESG (which is a directed graph) and a similarity subgraph (which is an undirected graph). They both share the same set of nodes, a node is defined as a spliced-exon (=s-exon), a s-exon is a set of aligned exonic sequences belonging to different species, in other words it is a MSA. In the ESG, there exist an edge between two nodes if there exist a transcript such that the sequences associated to each node co-occurs in that transcript. In the similarity subgraph, there exist an edge between two nodes if the p-value between those nodes are below some threshold (set to 0.001), in other words if the consensus sequence for each node are similar.
Thus an ASRU is a set of instances where an instance is composed of one or multiple s-exons.

A transcript is a path in the ESG. A translated transcript is a protein isoform thus the set of protein isoforms for a given gene is a subset of the set of paths of an ESG.
Hence an ESG represents the whole transcript variability observed in a set of species [1].

The notion of ASRU gives some highlights on the alternative usage of s-exons, hence on functional diversity generated by alternatively spliced events.


References
[1] : https://doi.org/10.1101/2020.11.14.382820
