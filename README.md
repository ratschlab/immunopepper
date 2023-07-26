# ImmunoPepper

**ImmunoPepper** is a software tool for the detection of neoantigens from a splicing graph. It generates the set of all theoretically feasible peptide sequences (or kmers) through direct translation of all walks along the graph.

The splice graph is generated using the SplAdder pipeline. SplAdder is a software tool for the detection of alternative splicing events from RNA-Seq data.
This splicing graph is generated from one or several alignment files derived from RNA-Seq data in bam format, as well as from a reference annotation file.
*Note: The bam file should be accompanied by a bam index file (.bai).*
For more information on how to incorporate this data to generate the input splice graph check out the [build mode of SplAdder](https://spladder.readthedocs.io/en/latest/spladder_modes.html)

The peptide set obtained by using this software tool can be personalized with germline and somatic variants, and takes all exon-exon junctions present in the splicing graph (even the ones not part of the reference annotation, but present in the given RNA-Seq sample) into account. The comprehensive set of peptides can be used subsequently for further downstream analyses, such as domain annotation or computational immunology.

For more information of the tool and its different modes check the online documentation: https://immunopepper.readthedocs.io/en/latest/
