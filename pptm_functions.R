

get_pptm <- function(x){
	x_total_count_sum <- sum(
		x$T0
		)
	T0_pptm <- x$T0 / (x_total_count_sum / 10^7)
	
	x_total_count_sum <- sum(
                x$Drug_1
                )
        Drug_1_pptm <- x$Drug_1 / (x_total_count_sum / 10^7)

	x_total_count_sum <- sum(
                x$Drug_2
                )
        Drug_2_pptm <- x$Drug_2 / (x_total_count_sum / 10^7)

	x_total_count_sum <- sum(
                x$DMSO_1
                )
        DMSO_1_pptm <- x$DMSO_1 / (x_total_count_sum / 10^7)

	x_total_count_sum <- sum(
                x$DMSO_2
                )
        DMSO_2_pptm <- x$DMSO_2 / (x_total_count_sum / 10^7)
	
	sgrna <- x$sgRNA
	gene <- x$Gene
	
	x_with_pptm <- cbind(
		sgrna,
		gene,
		T0_pptm,
                DMSO_1_pptm,
		DMSO_2_pptm,
		Drug_1_pptm,
		Drug_2_pptm
                )

	return(x_with_pptm)
}
