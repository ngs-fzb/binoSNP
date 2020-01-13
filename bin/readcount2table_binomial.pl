#!/usr/bin/perl
use Statistics::R;
use List::Util qw/sum/;

$num_args = $#ARGV + 1;

if ($num_args != 3) {
   print "\n Usage: readcount2table.pl /path/to/output/out.table /path/to/in.txt /path/to/Snps_RefuAlt.table\n";
    exit;
}
open(Fout,">$ARGV[0]") or die "\n\n\"Offnen von $ARGV[0] nicht m\"oglich\n\n\n";
open(Fin,"<$ARGV[1]") or die "\n\n\"Offnen von $ARGV[1] nicht m\"oglich\n\n\n";
open(Fin2,"<$ARGV[2]") or die "\n\n\"Offnen von $ARGV[2] nicht m\"oglich\n\n\n";

@ref=<Fin2>;#Zeilen der RefuAlt Tabelle in einem Array Speichern
print Fout "POS\tREF\tALT\tDP\t#A\tQual(A)\tFreq(A)\t#C\tQual(C)\tFreq(C)\t#G\tQual(G)\tFreq(G)\t#T\tQual(T)\tFreq(T)\tType\tFreq(ALT)\tQual(ALT)\tp-value\tAntibiotic\tComment\n";

while($line=<Fin>){
@array=split("\t", $line); #Jede Zeile des bam-readcount Datei nach Tabs auftrennen
  unless(defined $h{$array[1]}){ #Doppeleintr\"age nur einmal verarbeiten
@As=split(":", $array[5]); #Informationen zur Base A nach : auftrennen
@Cs=split(":", $array[6]);#Informationen zur Base C nach : auftrennen
@Gs=split(":", $array[7]);#Informationen zur Base G nach : auftrennen
@Ts=split(":", $array[8]);#Informationen zur Base T nach : auftrennen
@DEL=split(":", $array[10]);#Informationen ob es Deletion gibt nach : auftrennen
@index= grep(/^$array[1]/i,@ref); #Position der bam-eadcount-Datei mit der Position der RefuAlt-Tabelle matchen


$anzahl=@index; #Pr\"ufvariable, ob es mehrere alternative Basen und somit Doppeleintr\"age in der Resiliste gibt
if($array[3] == 0){#Abfangen der illegal division by zero
@zeile=split("\t", $index[0]);
		  chomp($zeile[3]);
print Fout "$zeile[0]\t$zeile[1]\t$zeile[2]\t$array[3]\t$As[1]\t$As[3]\t0\t$Cs[1]\t$Cs[3]\t0\t$Gs[1]\t$Gs[3]\t0\t$Ts[1]\t$Ts[3]\t-\t-\t-\t-\t-\t$zeile[3]\t\n";
}
elsif(defined $DEL[0]){#Abfrage, ob es sich um eine Deletion handelt
    @zeile=split("\t", $index[0]);
		  chomp($zeile[3]);
    print Fout "$zeile[0]\t$zeile[1]\t$zeile[2]\t$array[3]\t$As[1]\t$As[3]\t0\t$Cs[1]\t$Cs[3]\t0\t$Gs[1]\t$Gs[3]\t0\t$Ts[1]\t$Ts[3]\t0\tDEL\t-\t-\t-\t$zeile[3]\t$DEL[0] in $DEL[1] reads\n";}
else{
$freqA=$As[1]/$array[3];#Frequenz der Basen bestimmen
$freqC=$Cs[1]/$array[3];
$freqG=$Gs[1]/$array[3];
$freqT=$Ts[1]/$array[3];

my $R = Statistics::R->new();
$R->startR ;
#nun wird die Pr\"ufvariable abgefragt und je nach Anzahl der Eintr\"age in der Resiliste gibt es 1, 2 oder 3 Eintr\"age in der Tabelle
if ($anzahl == 1) {@zeile=split("\t", $index[0]);
		  chomp($zeile[3]);
		  if ($zeile[2] eq "A"){$freqALT=$freqA;
				      $qualALT=$As[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($As[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "C") {$freqALT=$freqC;
				      $qualALT=$Cs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Cs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "G") {$freqALT=$freqG;
				      $qualALT=$Gs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Gs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "T") {$freqALT=$freqT;
				      $qualALT=$Ts[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Ts[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		  else {print "Fehler Else1\n";}
	print Fout "$zeile[0]\t$zeile[1]\t$zeile[2]\t$array[3]\t$As[1]\t$As[3]\t$freqA\t$Cs[1]\t$Cs[3]\t$freqC\t$Gs[1]\t$Gs[3]\t$freqG\t$Ts[1]\t$Ts[3]\t$freqT\tSNP\t$freqALT\t$qualALT\t$Q\t$zeile[3]\t\n";
		    $h{$array[1]}=1;
}
  elsif ($anzahl == 2) {
		@zeile=split("\t", $index[0]);
		@zeile2=split("\t", $index[1]);
		chomp($zeile[3]);
		  if ($zeile[2] eq "A"){$freqALT=$freqA;
				      $qualALT=$As[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($As[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "C") {$freqALT=$freqC;
				      $qualALT=$Cs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Cs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "G") {$freqALT=$freqG;
				      $qualALT=$Gs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Gs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "T") {$freqALT=$freqT;
				      $qualALT=$Ts[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Ts[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		  else {print "Fehler Else2\n";}
		 chomp($zeile2[3]);
		  if ($zeile2[2] eq "A"){$freqALT2=$freqA;
				      $qualALT2=$As[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($As[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q2=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile2[2] eq "C") {$freqALT2=$freqC;
				      $qualALT2=$Cs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Cs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q2=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile2[2] eq "G") {$freqALT2=$freqG;
				      $qualALT2=$Gs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Gs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q2=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile2[2] eq "T") {$freqALT2=$freqT;
				      $qualALT2=$Ts[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Ts[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q2=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		  else {print "Fehler Else3\n";}
	print Fout "$zeile[0]\t$zeile[1]\t$zeile[2]\t$array[3]\t$As[1]\t$As[3]\t$freqA\t$Cs[1]\t$Cs[3]\t$freqC\t$Gs[1]\t$Gs[3]\t$freqG\t$Ts[1]\t$Ts[3]\t$freqT\tSNP\t$freqALT\t$qualALT\t$Q\t$zeile[3]\t\n";
	print Fout "$zeile2[0]\t$zeile2[1]\t$zeile2[2]\t$array[3]\t$As[1]\t$As[3]\t$freqA\t$Cs[1]\t$Cs[3]\t$freqC\t$Gs[1]\t$Gs[3]\t$freqG\t$Ts[1]\t$Ts[3]\t$freqT\tSNP\t$freqALT2\t$qualALT2\t$Q2\t$zeile2[3]\t\n";
	$h{$array[1]}=1;
   }
  elsif ($anzahl == 3) {
		@zeile=split("\t", $index[0]);
		@zeile2=split("\t", $index[1]);
		@zeile3=split("\t", $index[2]);
		chomp($zeile[3]);
		  if ($zeile[2] eq "A"){$freqALT=$freqA;
				      $qualALT=$As[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($As[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "C") {$freqALT=$freqC;
				      $qualALT=$Cs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Cs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "G") {$freqALT=$freqG;
				      $qualALT=$Gs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Gs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile[2] eq "T") {$freqALT=$freqT;
				      $qualALT=$Ts[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Ts[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    else {print "Fehler Else4\n";}
		chomp($zeile2[3]);
		  if ($zeile2[2] eq "A"){$freqALT2=$freqA;
				      $qualALT2=$As[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($As[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q2=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile2[2] eq "C") {$freqALT2=$freqC;
				      $qualALT2=$Cs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Cs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q2=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile2[2] eq "G") {$freqALT2=$freqG;
				      $qualALT2=$Gs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Gs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q2=$R->get('p');#mittlerer zu erwartender Sequenzierfehler
				      }
		    elsif ($zeile2[2] eq "T") {$freqALT2=$freqT;
				      $qualALT2=$Ts[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Ts[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q2=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    else {print "Fehler Else5\n";}
		 chomp($zeile3[3]);
		  if ($zeile3[2] eq "A"){$freqALT3=$freqA;
				      $qualALT3=$As[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($As[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q3=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile3[2] eq "C") {$freqALT3=$freqC;
				      $qualALT3=$Cs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Cs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q3=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile3[2] eq "G") {$freqALT3=$freqG;
				      $qualALT3=$Gs[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Gs[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q3=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		    elsif ($zeile3[2] eq "T") {$freqALT3=$freqT;
				      $qualALT3=$Ts[3];
				      $R->run(qq'quality<-c($As[3],$Cs[3],$Gs[3],$Ts[3])
				      meanerror<-10^(-mean(quality[which(quality>0)])/10)
				      p<-binom.test($Ts[1],$array[3],meanerror,alternative="greater", conf.level=0.95)');
				      $R->run(q'p<-p$p.value');
				      $Q3=$R->get('p');#mittlerer zu erwartender Fehler
				      }
		  else {print "Fehler Else6\n";}
	print Fout "$zeile[0]\t$zeile[1]\t$zeile[2]\t$array[3]\t$As[1]\t$As[3]\t$freqA\t$Cs[1]\t$Cs[3]\t$freqC\t$Gs[1]\t$Gs[3]\t$freqG\t$Ts[1]\t$Ts[3]\t$freqT\tSNP\t$freqALT\t$qualALT\t$Q\t$zeile[3]\t\n";
	print Fout "$zeile2[0]\t$zeile2[1]\t$zeile2[2]\t$array[3]\t$As[1]\t$As[3]\t$freqA\t$Cs[1]\t$Cs[3]\t$freqC\t$Gs[1]\t$Gs[3]\t$freqG\t$Ts[1]\t$Ts[3]\t$freqT\tSNP\t$freqALT2\t$qualALT2\t$Q2\t$zeile2[3]\t\n";
	print Fout "$zeile3[0]\t$zeile3[1]\t$zeile3[2]\t$array[3]\t$As[1]\t$As[3]\t$freqA\t$Cs[1]\t$Cs[3]\t$freqC\t$Gs[1]\t$Gs[3]\t$freqG\t$Ts[1]\t$Ts[3]\t$freqT\tSNP\t$freqALT3\t$qualALT3\t$Q3\t$zeile3[3]\t\n";
	$h{$array[1]}=1;
  }
else{print "Was ist mit Position $array[1]?\n"}  #Fehler abfangen

 $R->stopR() ;} } #unless
} #while
