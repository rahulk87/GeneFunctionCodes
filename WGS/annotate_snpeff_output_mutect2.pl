#!/usr/bin/perl
open(IN,$ARGV[0]);
	while($line=<IN>)
	{
	chomp $line;
	@sp=split('\t',$line);
		@gt=split(':',$sp[5]);
		if($gt[0] eq '0/0')
                        {
                        $zyf="Hom";
                        }
                        elsif($gt[0] eq '0/1')
                        {
                        $zyg="Het";
                        }
                        elsif($gt[0] eq '1/1')
                        {
                        $zyg="Hom";
                        }
                        elsif($gt[0] eq '1/2')
                        {
                        $zyg="Het";
                        }
		if($sp[17] eq '')
		{
		print "$sp[0]\t$sp[1]\t$sp[2]\t$sp[3]\t$sp[4]\t$sp[5]\t$sp[6]\t$sp[7]\t$sp[8]\t$sp[9]\t$sp[10]\t$sp[11]\t$sp[12]\t$sp[13]\t$sp[14]\t$sp[15]\t$sp[16]\tna\t$zyg\tNon-coding\tna\n";
		#print "$line\tNon-coding\tna\n";
		}
		elsif($sp[17] ne '')
		{
		$res=`grep -w '^$sp[14]' /scratch/DBC/GENFUNC/NGS_Projects/scripts/CGC_genes_191016.txt`;chomp $res;	
			if($res ne '')
         	        {
                	print "$line\t$zyg\tCoding\tCGC\n";
                	}
               		else
                	{
                	print "$line\t$zyg\tCoding\tNon-CGC\n";
                	}
		}
	}
