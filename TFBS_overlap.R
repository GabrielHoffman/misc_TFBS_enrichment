# Gabriel Hoffman

# Compute TFBS enrichment and p-values

# in bash
# ml openssl

library(GenomicRanges)
library(data.table)
library(regioneR)

# readEpiMap peaks
library(synapser)
library(genomation)
library(parallel)
synLogin()
epimapPeaks = readBed( synGet('syn8080422')$path )

# just keep chr10 for SPEED
epimapPeaks = epimapPeaks[seqnames(epimapPeaks) == "chr10"]

# read TFBS
dfTFBS = fread(cmd="zcat /sc/orga/projects/psychencode/data/external/jaspar2018/JASPAR2018_hg38_all_chr.bed.gz | awk '{if($5>600) print $0}' | head -n 100000")
grTFBS = with(dfTFBS, GRanges(V1, IRanges(V2, V3), names=V4, score=V5))

# just keep chr10 for SPEED
grTFBS = grTFBS[seqnames(grTFBS) == "chr10"]


# list of GRanges, one for each TF
#### Keep only 50 TFs for SPEED
grTFBS_list = lapply( unique( grTFBS$names )[1:50], function( tfbs ){
 grTFBS[grTFBS$names ==tfbs]
})
grTFBS_list = as(grTFBS_list, "GRangesList")
names(grTFBS_list) = lapply(grTFBS_list, function(x) x$names[1])

# for each TF, compute number of overlpas in epimapPeaks
count_TF_overlaps = sapply( grTFBS_list, function( gr ){
 numOverlaps( gr, epimapPeaks, count.once=TRUE )
})

# randomizations
# compute overlap of TFBS with random intervals
# need to customize the random peak set to the hypothesis we want to test
rnd_tf_matrix = mclapply( 1:500, function(x){

	rndPeaks = randomizeRegions(epimapPeaks)

	rnd_TF_overlaps = sapply( grTFBS_list, function( gr ){
	 	numOverlaps( gr, rndPeaks, count.once=TRUE )
	})
	rnd_TF_overlaps
}, mc.cores=12)
rnd_tf_matrix = do.call("rbind", rnd_tf_matrix)

# compute statistics for each TF
pseudocount = 0.25
result = lapply( names(count_TF_overlaps), function(tf){

	# get p-value based on permuations
	p.perm = sum(rnd_tf_matrix[,tf] > count_TF_overlaps[tf]) / nrow(rnd_tf_matrix) 

	# compute enrichment
	tab = c(count_TF_overlaps[tf], 
			length(epimapPeaks), 
			median(rnd_tf_matrix[,tf]),		
			length(epimapPeaks)	)
	tab = matrix(tab,2,2)
	fit = fisher.test(tab)

	# pnorm( count_TF_overlaps[tf], mean(rnd_tf_matrix[,tf]), sd=sd(rnd_tf_matrix[,tf]), lower.tail=FALSE )
	zscore = (count_TF_overlaps[tf] - mean(rnd_tf_matrix[,tf])) / sd(rnd_tf_matrix[,tf])
	
	data.frame( TF=tf, 
				log2OR = log2(fit$estimate),
				ci.lower = log2(fit$conf.int[1]),
				ci.upper = log2(fit$conf.int[2]),
				p.value = fit$p.value,
				p.perm,
				zscore = zscore,
				p.zscore = pnorm(zscore, lower.tail=FALSE))
})
result = do.call("rbind", result)
rownames(result) = c()

result

library(ggplot2)

# plot log2 OR
ggplot(result, aes(TF,log2OR)) + geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), color="grey40") + geom_point() + coord_flip() + theme_bw() + geom_hline(yintercept=0) 


# plot zscores
# ggplot(result, aes(TF,zscore)) + geom_point() + coord_flip() + theme_bw() + geom_vline(xintercept=0)







