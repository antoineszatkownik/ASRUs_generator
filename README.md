# ASRUs_generator

ASRU stands for Alternatively Spliced Repetitive Unit. An ASRU is a connected component in the similarity subgraph of an Evolutionnary Splicing Graph (=ESG). 
For a given gene we associate to it two graphs, an ESG and a similarity subgraph. They both share the same set of nodes, a node is defined by a spliced-exon (=s-exon), a s-exon is a set of aligned exonic sequences belonging to different species, in other words it is a MSA. In the ESG, there exist an edge between two nodes if there exist a transcript such that the sequences associated to each node co-occurs in that transcript. In the similarity subgraph, there exist an edge between two nodes if the p-value between those nodes are below some threshold (set to 0.001), in other words if the consensus sequence for each node are similar.
