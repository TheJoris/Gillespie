#!/usr/bin/perl
#Give project name of KaiModel Gillespie output as only argument.
# NOT including '.0.' part!


# OUTPUT is:

# <project_name>.CA       : KaiC seq. with KaiA. (On CII)
# <project_name>.CAB      : KaiC seq. with KaiB seq. KaiA.
# <project_name>.CB       : KaiC seq. with KaiB.
# <project_name>.Cfree    : Non sequestered KaiC (norm + inactive, int over p).
# <project_name>.tmin.Np  : -
# <project_name>.Puw      : Fraction of hexamers with phosho. level > 0.
# <project_name>.Pw       : Av. phospho level of hexamers: 1/(6 * CT) sum_i (p * C_p_i).

$|=1;


print "argv is @ARGV\n";
$Name    = @ARGV[0];

$twrite = 1900;

$nf = 2;
$file[0] = "$Name.0.A";
$file[1] = "$Name.0.B";

for ($f= 0; $f<$nf; $f++) {
    
    print "We read in $file[$f]\n";
    
    open(FILE,"$file[$f]");
    $c = 0; $t = 0;
    while(<FILE>) {
	      @data = split(' ',$_);
	      if ($c>0) {
	          $T[$f][$t] = @data[0];
	          if ($f==0) {
		            $A[$t] = @data[1];
		            if ($t == $twrite) {
		                open(FILE1,">$file[$f].$t");
		                printf FILE1 "%f %d\n",$T[0][$t],$A[$t];
		                close(FILE1);
		            }
	          }
	          elsif ($f==1) {
		            $B[$t] = @data[1];
		            if ($t == $twrite) {
		                open(FILE1,">$file[$f].$t");
		                printf FILE1 "%f %d\n",$T[0][$t],$B[$t];
		                close(FILE1);
		            }
	          }
	          elsif ($f==2) {
		            $AB[$t] = @data[1];
		            if ($t == $twrite) {
		                open(FILE1,">$file[$f].$t");
		                printf FILE1 "%f %d\n",$T[0][$t],$AB[$t];
		                close(FILE1);
		            }
	          }
	          $t++;
	      }
	      $c++;
    }
    close(FILE);
}

$tt = $t;
for ($t = 0; $t < $tt; $t++) {
    $AT[$t] = $A[$t]; # + $AB[$t];
    $BT[$t] = $B[$t]; # + $AB[$t];
}

for ($p = 0; $p<7; $p++) {

    $nf = 8;    
    $file[0] = "$Name.0.C$p";
    $file[1] = "$Name.0.C$p"."A";
    $file[2] = "$Name.0.C$p"."i";
    $file[3] = "$Name.0.C$p"."iB";
    $file[4] = "$Name.0.C$p"."iB2";
    $file[5] = "$Name.0.C$p"."iAB";
    $file[6] = "$Name.0.C$p"."iAB2";
    $file[7] = "$Name.0.C$p"."iA2B2";

    for ($f = 0; $f < $nf; $f++) {

	      print "We read in $file[$f]\n";

	      open(FILE,"$file[$f]");
	      $c = 0; $t = 0;
	      while(<FILE>) {
	          @data = split(' ',$_);
	          if ($c>0) {
		            $T[$f][$t] = @data[0];
		            $X[$f][$t] = @data[1];
		            if ($t == $twrite) {
		                open(FILE1,">$file[$f].$t");
		                printf FILE1 "%f %d\n",$T[0][$t],$X[$f][$t];
		                close(FILE1);
		            }
		            $t++;
	          }
	          $c++;
	      }
	      close(FILE);
	      close(FILE1);
    }
    print "T is $t\n";

    $tt = $t;
    for ($t=0; $t<$tt; $t++) {
	      $C[$p][$t] = 0;
	      $AT[$t] += $X[1][$t] + $X[5][$t] + $X[6][$t] + 2 * $X[7][$t];
	      $BT[$t] += $X[3][$t] + $X[5][$t] + 2*($X[4][$t]+$X[6][$t]+$X[7][$t]);
	      for ($f = 0; $f < $nf; $f++) {
	          $C[$p][$t] += $X[$f][$t];
	      }
	      $Cf[$p][$t]  = $X[0][$t] + $X[2][$t];
	      $CA[$p][$t]  = $X[1][$t];
	      $CB[$p][$t]  = $X[3][$t] + $X[4][$t];
	      $CAB[$p][$t] = $X[5][$t] + $X[6][$t] + $X[7][$t];
  }
}

$Pwmin = 1.;
for ($t = 0; $t < $tt; $t++) {

    $Pw[$t] = 0; $Puw[$t] = 0; 
    $CT[$t] = $C[0][$t]; 
    $CfT[$t] = $Cf[0][$t];
    $CAT[$t] = $CA[0][$t];
    $CBT[$t] = $CB[0][$t];
    $CABT[$t] = $CAB[0][$t];
    for ($p=1; $p < 7; $p++) {
        $CT[$t]  += $C[$p][$t];
        $CfT[$t] += $Cf[$p][$t];
        $CAT[$t] += $CA[$p][$t];
        $CBT[$t] += $CB[$p][$t];
        $CABT[$t] += $CAB[$p][$t];
        $Ctest = $CfT[$t] + $CAT[$t] + $CBT[$t] + $CABT[$t];
        $Pw[$t]  += $p * $C[$p][$t];
        $Puw[$t] += $C[$p][$t];
    }
    print "$t: CT is $CT[$t] or $Ctest; AT is $AT[$t]; BT is $BT[$t]\n";
		if ($CT[$t] > 0){
	    $Puw[$t] /= $CT[$t];
  	  $Pw [$t] /= (6 * $CT[$t]);
		}
		else
		{
	    $Puw[$t] = 0;
  	  $Pw [$t] = 0;
		}

    if ($Pw[$t] < $Pwmin && $T[0][$t] > 50.) {
  	    $tmin = $t;
	      $Pwmin = $Pw[$t];
    }
}

print "Pwmin is $Pwmin\n";
print "tmin is $tmin\n";
print "Tmin is $T[0][$tmin]\n";

$file = "$Name.$tmin.Np\n";
open(FILE,">$file");
printf FILE "phospho - # of KaiC in phospho. states [0,6] \n";
for ($p = 0; $p < 7; $p++) {
    printf FILE "%d\t%d\t0\n",$p,$C[$p][$tmin];
}
close(FILE);

$file = "$Name.Pw";
open(FILE,">$file");
printf FILE "time - sum_i i * C_i /(6 * CT) (i>0) \n";
for ($t=0;$t<$tt;$t++) {
    printf FILE "%f\t%f\t0\n",$T[$0][$t],$Pw[$t];
}
close(FILE);

$file = "$Name.Puw";
open(FILE,">$file");
printf FILE "time - sum_i C_i /CT (i>0) \n";
for ($t=0;$t<$tt;$t++) {
    printf FILE "%f\t%f\t0\n",$T[$0][$t],$Puw[$t];
}
close(FILE);

$file = "$Name.Cfree";
open(FILE,">$file");
printf FILE "time - non seq. KaiC \n";
for ($t=0;$t<$tt;$t++) {
    printf FILE "%f\t%f\t0\n",$T[$0][$t],$CfT[$t];
}
close(FILE);

$file = "$Name.CA";
open(FILE,">$file");
printf FILE "time - KaiC-KaiA (on CII) \n";
for ($t=0;$t<$tt;$t++) {
    printf FILE "%f\t%f\t0\n",$T[$0][$t],$CAT[$t];
}
close(FILE);

$file = "$Name.CB";
open(FILE,">$file");
printf FILE "time - KaiC-KaiB \n";
for ($t=0;$t<$tt;$t++) {
    printf FILE "%f\t%f\t0\n",$T[$0][$t],$CBT[$t];
}
close(FILE);

$file = "$Name.CAB";
open(FILE,">$file");
printf FILE "time - KaiC-KaiB-KaiA \n";
for ($t=0;$t<$tt;$t++) {
    printf FILE "%f\t%f\t0\n",$T[$0][$t],$CABT[$t];
}

$file = "$Name.CT";
open(FILE,">$file");
printf FILE "time - Total KaiC \n";
for ($t=0;$t<$tt;$t++) {
    printf FILE "%f\t%f\t0\n",$T[$0][$t],$CT[$t];
}

close(FILE);


