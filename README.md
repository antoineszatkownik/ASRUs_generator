# Context

ASRU stands for Alternatively Spliced Repetitive Unit. An ASRU is a connected component in the similarity subgraph of an Evolutionnary Splicing Graph (=ESG, [1]). 
For a given gene we associate to it two graphs, an ESG (which is a directed graph) and a similarity subgraph (which is an undirected graph). They both share the same set of nodes, a node is defined as a spliced-exon (=s-exon), a s-exon is a set of aligned exonic sequences belonging to different species, in other words it is a MSA. In the ESG, there exist an edge between two nodes if there exist a transcript such that the sequences associated to each node co-occurs in that transcript. In the similarity subgraph, there exist an edge between two nodes if the p-value between those nodes are below some threshold (set to 0.001), in other words if the consensus sequence for each node are similar.
Thus an ASRU is a set of instances where an instance is composed of one or multiple s-exons.

A transcript is a path in the ESG. A translated transcript is a protein isoform thus the set of protein isoforms for a given gene is a subset of the set of paths of an ESG.
Hence an ESG represents the whole transcript variability observed in a set of species [1].

The notion of ASRU gives some highlights on the alternative usage of s-exons, hence on functional diversity generated by alternatively spliced events.

# About the code

For a given gene, or a list of genes it can generates two .csv file, one about the ASRUs and the other one about the instances of the ASRUS, they are of the given form : <br />

gene,uniteRepetee,instances,max,min,moy,median,ecartType,evenements <br />
SLC6A12,"{'1_1', '1_7'}",2,68,37,52.5,52.5,15.5,[3] <br />
SLC6A12,"{'5_2', '5_0/5_1'}",2,48,47,47.5,47.5,0.5,[2]

instance,taille,nombre,UniteRepetee,gene <br />
1_1,68,1,"{'1_1', '1_7'}",SLC6A12 <br />
1_7,37,1,"{'1_1', '1_7'}",SLC6A12 <br />
5_2,47,1,"{'5_2', '5_0/5_1'}",SLC6A12 <br />
5_0/5_1,48,2,"{'5_2', '5_0/5_1'}",SLC6A12

where max is the maximum in length, and the events column are those events (in fact the rank of the event, from highly conserved to less conserved) where the s-exon(s) played a role in.




# How to use

The data used for this code are provided in the "data" folder. Prior to running the code you should not forget to rename the paths.
Then all you have to do is remove the hash sign at either of the following line (present at the end of the code) : <br />
<br />
#writecsv_ASRU()       
#writecsv_instances() <br />

writecsv_ASRU() will write the csv about the ASRUs, writecsv_instances() will write the csv about the instances.

***Note:***
The variable "earth" is the list of all gene's name highlighted by ThorAxe [1], it contains 2190 names of gene. If you are interested in a specific set of genes, then replace it by a list of the names of those genes.

***References*** <br />
[1] : https://doi.org/10.1101/2020.11.14.382820
