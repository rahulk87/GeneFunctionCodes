library(ggplot2)
library(ggdendro)
library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.1000genomes.hs37d5)

# mutect_output.txt is joined mutect file
sca_vr = readMutect("mutect_output.txt")
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.1000genomes.hs37d5)
sca_mm = motifMatrix(sca_motifs, group = "sampleNames", normalize = TRUE)
png("somatic_sig.png",height=3000, width=3000, res=400)
plotMutationSpectrum(sca_motifs, "sampleNames")
dev.off()

########
#  Number of signatures
n_sigs = 2:8
gof_nmf = assessNumberSignatures(sca_mm, n_sigs, nReplicates = 5)
png("nmf.png",height=3000, width=3000, res=400)
plotNumberSignatures(gof_nmf)
dev.off()

#######
# Plot signature
n_sigs = 5
sigs_nmf = identifySignatures(sca_mm, n_sigs, nmfDecomposition)
png("signature_bar.png",height=3000, width=3000, res=400)
plotSignatures(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")
dev.off()

#######
# Plot signature
png("signature_heat.png",height=3000, width=3000, res=400)
plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")
dev.off()

#######
# Plot observed spectrum
png("signature_obs.png",height=3000, width=3000, res=400)
plotObservedSpectrum(sigs_nmf)
dev.off()

#######
# Plot fitted spectrum
png("signature_fit.png",height=3000, width=3000, res=400)
plotFittedSpectrum(sigs_nmf)
dev.off()


#######
# Plot signature contribution in each sample
png("signature_contribution_heat.png",height=3000, width=3000, res=400)
plotSampleMap(sigs_nmf)
dev.off()


#######
# Plot signature contribution in each sample (bar plot)
png("signature_contribution_bar.png",height=3000, width=3000, res=400)
plotSamples(sigs_nmf)
dev.off()



#######
# Plot clustering of samples
clu_motif = clusterSpectrum(sca_mm, "motif")
png("cluster.png",height=3000, width=3000, res=400)
ggdendrogram(clu_motif, rotate = TRUE)
dev.off()
