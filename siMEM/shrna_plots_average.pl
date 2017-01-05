# This perl script generates the average dropout plots from COLT2 data. User has to provide the cell line classification file (e.g. RB1_Null or RB1_WT here).

#!/usr/bin/perl
open(IN5,"../data/version9_classification.txt");
while($line5=<IN5>)
{
chomp $line5;
@sp5=split('\t',$line5);
`grep -w '$sp5[0]' ../data/cell_line_position_in_shrna_data.txt |cut -f1 >>pos`;
}
`perl -pi -e 's/\n/,/g' pos`;
`perl -pi -e 's/,\$//g' pos`;
open(IN6,"pos");
$line6=<IN6>;chomp $line6;
#print "$line6\n";
`cut -f1,2,$line6 ../data/log2_matrix_with_gene_symbol_mod.txt >log2_matrix`;
`head -1 log2_matrix >head`;

open(IN,"../data/list");
	while($line=<IN>)
	{
	chomp $line;print "Running for $line....\n";
	open(FA,">>plot_input_$line");
	print FA "Cell_line\tTime_points\tshRNA_Dropout\tRB1_Status\n";
	`head -2 log2_matrix >tmp`;
	`grep -w '$line' log2_matrix >>tmp`;
	`R CMD BATCH ../codes/transpose.R`;
	`grep -v '""' y |perl -pi -e 's/"//g' |perl -pi -e 's/,/\t/g' |grep -v 'symbol' |grep -v 'trcn_id' |perl -pi -e 's/breast_//g' >tmp1`; 
	#`grep -v '""' y |perl -pi -e 's/"//g' |perl -pi -e 's/,/\t/g' |grep -v 'symbol' |grep -v 'trcn_id' |perl -pi -e 's/breast_//g' |perl -pi -e 's/_t0_a//g' |perl -pi -e 's/_t0_b//g' |perl -pi -e 's/_t0_c//g' |perl -pi -e 's/_t1_a//g' |perl -pi -e 's/_t1_b//g' |perl -pi -e 's/_t1_c//g' |perl -pi -e 's/_t2_a//g' |perl -pi -e 's/_t2_b//g' |perl -pi -e 's/_t2_c//g' >tmp1`;
	#`grep -w 'tnbc' tmp1 |grep -v 'jimt1' |grep -v 'mcf10a' >tmp2`; ## RB1null_mut_CNV
	#`grep -w 'tnbc' tmp1 >tmp2`;
	`cut -f1 tmp1 |cut -d '_' -f1 >tmp3`;
	`cut -f1 tmp1 |cut -d '_' -f1 |sort -u >cells`;
	`paste tmp3 tmp1 >tmp4`;
		open(FC,">>for_ttest");
		print FC "cells\tstatus\tt0\tt1\tt2\n";
		open(IN1,"cells");
		while($line1=<IN1>)
		{
		chomp $line1;
		$status=`grep -w $line1 ../data/version9_classification.txt |cut -f2`;chomp $status;
		$t0=$t1=$t2=0;
		### T0 ###
		`grep -w '$line1' tmp4 |grep 't0' |cut -f4- |perl -pi -e 's/\t/\n/g' >tmp5`;
		$c1=`wc -l tmp5 |perl -p -e 's/^ \\s+//g'|cut -d ' ' -f1`;chomp $c;
			open(IN2,"tmp5");
                        while($line2=<IN2>)
                        {
                        chomp $line2;
                        $t0=$t0+$line2;
                        }
		$avg_t0=$t0/$c1;
		print FA "$line1\tT0\t$avg_t0\t$status\n";
		### T1 ###
		`grep -w '$line1' tmp4 |grep 't1' |cut -f4- |perl -pi -e 's/\t/\n/g' >tmp5`;
                $c2=`wc -l tmp5 |perl -p -e 's/^ \\s+//g'|cut -d ' ' -f1`;chomp $c;
                        open(IN3,"tmp5");
                        while($line3=<IN3>)
                        {
                        chomp $line3;
                        $t1=$t1+$line3;
                        }
                $avg_t1=$t1/$c2;
                print FA "$line1\tT1\t$avg_t1\t$status\n";
		### T2 ###
                `grep -w '$line1' tmp4 |grep 't2' |cut -f4- |perl -pi -e 's/\t/\n/g' >tmp5`;
                $c3=`wc -l tmp5 |perl -p -e 's/^ \\s+//g'|cut -d ' ' -f1`;chomp $c;
                        open(IN4,"tmp5");
                        while($line4=<IN4>)
                        {
                        chomp $line4;
                        $t2=$t2+$line4;
                        }
                $avg_t2=$t2/$c3;
                print FA "$line1\tT2\t$avg_t2\t$status\n";
		print FC "$line1\t$status\t$avg_t0\t$avg_t1\t$avg_t2\n";
		}
	open(FB, ">>box.R");
	print FB "library(ggplot2)\n";
	print FB "pdf(\"$line.pdf\")\n";
	print FB "x<-read.csv(\"plot_input_$line\",sep=\"\\t\",header=TRUE)\n";
	print FB "ggplot(x, aes(x=Time_points, y=shRNA_Dropout, color=RB1_Status))+geom_boxplot(position=position_dodge(0.8))+geom_jitter(position=position_dodge(0.8)) + ggtitle(\"$line\") +  theme(plot.title = element_text(lineheight=.8, face=\"bold\"))\n";
	print FB "dev.off()\n";
	close FB;
	close FA;
	close FC;
	`R CMD BATCH box.R`;
	`rm tmp* transpose.Rout cells box.R*`;
	## T test ##
	open(FD,">>pvalue_$line");
	print FD "t0\tt1\tt2\n";
	for($i=3;$i<=5;$i++)
	{
	`grep 'WT' for_ttest |cut -f$i >wt`;
	`grep 'Null' for_ttest |cut -f$i >null`;
	`paste wt null >profile`;
	`Rscript ../../../ttest.R profile`;
	$p=`grep '"p.value"' y |cut -d ',' -f2`;chomp $p;
	print FD "$p\t";
	`rm wt null profile`;
	}
	`rm for_ttest`;
	`rm plot_input_* pvalue_*`;
	}
`rm head pos log2_matrix y`;
