#!/usr/bin/perl

print "# Generationn\tmostfit\tmaxfitdiv\tminfitdiv\tavgfitdiv\tvarfitdiv\tmaxphendiv\tminphendiv\tavgphendiv\tvarphendiv\tmaxhamdiv\tminhamdiv\tavghamdiv\tvarhamdiv\tTInowght\n";

$cycles=0;

while(<>) {
    chomp();

    if (/^STARTING/){
	$cycles++;
	$dataline = 0;
    }
    elsif (/^Starting/){
	@fields = split(/ /);
	$relevance = $fields[2];
    }
    elsif(length > 3) {
	@fields = split / /;
	$dataline++;
   	$realgen[$dataline]  = $fields[0];      # WARNING: FULL TEST NOT ON REGULAR GENERATION BASIS
	$mostfit[$dataline] += $fields[1];
	$maxfitdiv[$dataline]+= $fields[2];
	$minfitdiv[$dataline]+= $fields[3];
	$avgfitdiv[$dataline]+= $fields[4];
	$varfitdiv[$dataline]+= $fields[5];
	$maxphendiv[$dataline]+= $fields[6];
	$minphendiv[$dataline]+= $fields[7];
	$avgphendiv[$dataline]+= $fields[8];
	$varphendiv[$dataline]+= $fields[9];
	$maxhamdiv[$dataline]+= $fields[10];
	$minhamdiv[$dataline]+= $fields[11];
	$avghamdiv[$dataline]+= $fields[12];
	$varhamdiv[$dataline]+= $fields[13];
	$tinowght[$dataline]+= $fields[14];
    }
    elsif (length($_) < 3) {
	## compute averages for this batch:
	print "# data for relevance = $relevance\n";

	for ($i=1 ; $i < @mostfit; $i++){
	    $avg_mostfit = sprintf "%1.3f", ($mostfit[$i]/$cycles);
	    $avg_maxfitdiv = sprintf "%1.3f", ($maxfitdiv[$i]/$cycles); 
	    $avg_minfitdiv = sprintf "%1.3f", ($minfitdiv[$i]/$cycles); 
	    $avg_avgfitdiv = sprintf "%1.3f", ($avgfitdiv[$i]/$cycles); 
	    $avg_varfitdiv = sprintf "%1.3f", ($varfitdiv[$i]/$cycles); 
	    $avg_maxphendiv = sprintf "%1.3f", ($maxphendiv[$i]/$cycles); 
	    $avg_minphendiv = sprintf "%1.3f", ($minphendiv[$i]/$cycles); 
	    $avg_avgphendiv = sprintf "%1.3f", ($avgphendiv[$i]/$cycles); 
	    $avg_varphendiv = sprintf "%1.3f", ($varphendiv[$i]/$cycles); 
	    $avg_maxhamdiv = sprintf "%1.3f", ($maxhamdiv[$i]/$cycles); 
	    $avg_minhamdiv = sprintf "%1.3f", ($minhamdiv[$i]/$cycles); 
	    $avg_avghamdiv = sprintf "%1.3f", ($avghamdiv[$i]/$cycles); 
	    $avg_varhamdiv = sprintf "%1.3f", ($varhamdiv[$i]/$cycles); 
	    $avg_tinowght = sprintf "%1.3f", ($tinowght[$i]/$cycles);
	    #$dataline=$i*100;

	    print "$realgen[$i]\t$avg_mostfit\t$avg_maxfitdiv\t$avg_minfitdiv\t$avg_avgfitdiv\t$avg_varfitdiv\t$avg_maxphendiv\t$avg_minphendiv\t$avg_avgphendiv\t$avg_varphendiv\t$avg_maxhamdiv\t$avg_minhamdiv\t$avg_avghamdiv\t$avg_varhamdiv\t$avg_tinowght\n";
	}
	print "\n";
	# reset collectors
	for ($i=1 ; $i < @mostfit; $i++){
	    $mostfit[$i] = 0;
	    $maxfitdiv[$i] = 0;
	    $minfitdiv[$i] = 0;
	    $avgfitdiv[$i] = 0;
	    $varfitdiv[$i] = 0;
	    $maxphendiv[$i] = 0;
	    $minphendiv[$i] = 0;
	    $avgphendiv[$i] = 0;
	    $varphendiv[$i] = 0;
	    $maxhamdiv[$i] = 0;
	    $minhamdiv[$i] = 0;
	    $avghamdiv[$i] = 0;
	    $varhamdiv[$i] = 0;
	    $tinowght[$i] = 0;
	}
	$cycles=0;
    }

}

print "\n";
