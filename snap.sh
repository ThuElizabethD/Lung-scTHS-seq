#hg38_csize="/media/Scratch_SSD_Voyager/dinh/Genomes/bsos_hg38_bwa_index/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai"
hg38_size="hg38_chromsizes_only_xsms_notY.txt"
#sort bam file

#create snap file
snaptools snap-pre  \
--input-file=hALL.sorted.bam  \
--output-snap=hALL.merged.snap  \
--genome-name=hg38  \
--genome-size=$hg38_size  \
--min-mapq=30  \
--min-flen=0  \
--max-flen=1000  \
--keep-chrm=TRUE  \
--keep-single=TRUE  \
--keep-secondary=False  \
--overwrite=True  \
--min-cov=100  \
--verbose=True

# create bins x cells matrix /media/Home_Raid1/eli.duong/.local/bin 
snaptools snap-add-bmat	\
	--snap-file=hALL.merged.snap \
	--bin-size-lis 5000	\
	--verbose=True

