# circuclust

Circular alignment and clustering of nt sequences.

all-vs-all circular global alignment of query-vs-database:

<pre>
  circuclust -align query.fa \
	-db db.fa \
	-output alignments.txt \
	-tsvout hits.tsv \
</pre>

Cluster using UCLUST algorithm with circular global alignment,
  both plus and minus strands are checked,
  minimum fractional identity is specified by -id option:

<pre>
  circuclust -cluster seqs.fa \
 	-id 0.9 \
	-fastaout centroids.fa \
	-tsvout hits.tsv
</pre>

### Download pre-compiled binaries
Go to **Releases** at the github project home [https://github.com/rcedgar/circuclust/releases](https://github.com/rcedgar/circuclust/releases "https://github.com/rcedgar/circuclust/releases").

### Compile from source
Linux: The gcc compiler and make are required. Clone / download the repository, go to the `src/` sub-directory and type `make`. If successful, the binary will be in the `bin/` sub-directory.

Windows / Visual Studio: Load and build the `circuclust.sln` solution file.

Mac and Windows Cygwin: The Linux Makefile should work more or less, not tested so some tweaking may be needed.
