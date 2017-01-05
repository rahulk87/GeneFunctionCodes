# This perl script generates the hairpin wise dropout plots from COLT2 data. User has to provide the cell line classification file (e.g. RB1_Null or RB1_WT here).

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
	print FA "Cell_line\tTime_points\tshRNA_Dropout\tHairpin_ID\tRB1_Status\n";
	`head -2 log2_matrix >tmp`;
	`grep -w '$line' log2_matrix >>tmp`;
	`R CMD BATCH ../codes/transpose.R`;
	`grep -v '""' y |perl -pi -e 's/"//g' |perl -pi -e 's/,/\t/g' |grep -v 'symbol' |perl -pi -e 's/breast_//g' >tmp1`; 
	`cut -f1 tmp1 |cut -d '_' -f1 >tmp3`;
	`cut -f1 tmp1 |cut -d '_' -f1 |sort -u |grep -v 'trcn' >cells`;
	`paste tmp3 tmp1 >tmp4`;
	$hp=`head -1 tmp4`;chomp $hp;@hp1=split('\t',$hp);$hp2=scalar(@hp1);#print "$hp2\n";
		open(IN1,"cells");
		while($line1=<IN1>)
		{
		chomp $line1;#print "$line1....\n";
		$status=`grep -w $line1 ../data/version9_classification.txt |cut -f2`;chomp $status;
			for($i=4;$i<=$hp2;$i++)
			{
			$a=$i-1;#print"$i\n";
			$hp_id=@hp1[$a];chomp $hp_id;
			$t0=$t1=$t2=0;
			### T0 ###
			`grep -w '$line1' tmp4 |grep 't0' |cut -f$i >tmp5`;
			$c1=`wc -l tmp5 |perl -p -e 's/^ \\s+//g'|cut -d ' ' -f1`;chomp $c1;
				open(IN2,"tmp5");
                        	while($line2=<IN2>)
                        	{
                        	chomp $line2;
                        	$t0=$t0+$line2;
                        	}
                		$avg_t0=$t0/$c1;
			print FA "$line1\tT0\t$avg_t0\t$hp_id\t$status\n";
			### T1 ###
                        `grep -w '$line1' tmp4 |grep 't1' |cut -f$i >tmp5`;
                        $c2=`wc -l tmp5 |perl -p -e 's/^ \\s+//g'|cut -d ' ' -f1`;chomp $c1;
                                open(IN3,"tmp5");
                                while($line3=<IN3>)
                                {
                                chomp $line3;
                                $t1=$t1+$line3;
                                }
                                $avg_t1=$t1/$c2;
                        print FA "$line1\tT1\t$avg_t1\t$hp_id\t$status\n";
			### T2 ###
                        `grep -w '$line1' tmp4 |grep 't2' |cut -f$i >tmp5`;
                        $c3=`wc -l tmp5 |perl -p -e 's/^ \\s+//g'|cut -d ' ' -f1`;chomp $c1;
                                open(IN4,"tmp5");
                                while($line4=<IN4>)
                                {
                                chomp $line4;
                                $t2=$t2+$line4;
                                }
                                $avg_t2=$t2/$c3;
                        print FA "$line1\tT2\t$avg_t2\t$hp_id\t$status\n";
			}		
		}
	close FA;
	open(FB, ">>box.R");
        print FB "library(ggplot2)\n";
        print FB "pdf(\"$line.pdf\")\n";
        print FB "x<-read.csv(\"plot_input_$line\",sep=\"\t\",header=TRUE)\n";
        print FB "ggplot(x, aes(x = Time_points, y = shRNA_Dropout, color = RB1_Status)) + geom_boxplot() + facet_wrap(~ Hairpin_ID) + ggtitle(\"$line\") + geom_jitter(position=position_dodge(0.8)) + geom_boxplot(position=position_dodge(0.8)) + theme(plot.title = element_text(lineheight=.8, face=\"bold\"))\n";
	print FB "dev.off()\n";
        close FB;
        `R CMD BATCH box.R`;
	`rm tmp* cells transpose.Rout box.R* y plot_input_*`;
	}
`rm head pos log2_matrix`;
