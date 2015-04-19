<?php
$sGenotypeFile  = "";
$sExportPath = "";
$sMarkerSubsetFile = "";
$nPop1=0;
$nPop2=0;
$nPop3=0;

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-s':
            $sGenotypeFile   = array_shift($argv);
            break;
        case '-o':
            $sExportPath  = array_shift($argv);
            break;
		case '-i':
			$sMarkerSubsetFile  = array_shift($argv);
			break;
		case '-pop1':
			$nPop1  = array_shift($argv);
			break;
		case '-pop2':
			$nPop2  = array_shift($argv);
			break;
		case '-pop3':
			$nPop3  = array_shift($argv);
			break;
    }
}
 
/* print them */
echo("This is php. $sGenotypeFile, $sExportPath, $sMarkerSubsetFile.\n");

//read marker subset file:
$fMarkerSubsetFile = fopen($sMarkerSubsetFile, "r");
$arrMarkersToSample = array();

while(($sLn = fgets( $fMarkerSubsetFile   ))!==false) {
	//echo($sLn);
	$arrCols = explode("\t", $sLn);
	$nChr = $arrCols[0];
	$nIndex = $arrCols[1];
	if (!isset($arrMarkersToSample[$nChr])) {
		$arrMarkersToSample[$nChr] = array();
	}
	
	array_push( $arrMarkersToSample[$nChr] , intval($nIndex));
	
}
fclose($fMarkerSubsetFile );

$sPriorAlleleFreqFile = "$sExportPath/priorallelefreqs.txt";
$sExGenotypesFile = "$sExportPath/genotypes.txt";

 $fGenotypeFile = fopen($sGenotypeFile , "r"); // source
 
 $fPriorAlleleFreqFile = fopen($sPriorAlleleFreqFile , "w");
 $fExGenotypesFile = fopen($sExGenotypesFile , "w");
 
 $sLn = "";
 $bHeaderWritten = false;
 
 if ($fGenotypeFile == false ) {
	echo("Input file not found\n");
 }
 

 $arrAlleleCount_bigA_Pop1 = array();
 $arrAlleleCount_bigA_Pop2 =  array();
 $arrAlleleCount_littleA_Pop1 =  array();
 $arrAlleleCount_littleA_Pop2 =  array();
 $arrAlleleName =  array();
 $arrAlleleTemp = array();
 $nChrSide=0;
 

 
 //print(count($arrMarkersToSample[1]));
 //die();
// die();
$nCurrPop1=0;
$nCurrPop2=0;
$nCurrPop3=0;

 $nLnCount = -1;
 $bFirstPop3 = true;

 while(($sLn = fgets( $fGenotypeFile  ))!==false) {
	 $nLnCount++;
	$arrColumns = explode("\t", $sLn, 5);
	//echo("lINE:".$nLnCount.PHP_EOL);
	
	$nIndId = $arrColumns[0];
	$nPopId = $arrColumns[1];
	$nSex	= $arrColumns[2];
	
	
	//echo($arrColumns[0].PHP_EOL);
	$nMarkersToSampleCount = 0;
	
	if ($nPopId == 3 && $nLnCount ==0) {
		die("ERROR: no samples are present for the parental populations. This is required for exporting to admixmap. Rerun your simulations with parental populations!");
	}
	
	if ($nPopId  != 3) { // parental populations, used for calculating prior allele freq.
		if ($nPopId == 1) {
			$nCurrPop1++;
			if ($nCurrPop1>$nPop1) {
				continue;
			}
		}
		if ($nPopId == 2) {
			$nCurrPop2++;
			if ($nCurrPop2>$nPop2) {
				continue;
			}
		}
		
		$arrColumns = array_slice(explode("\t", $sLn) , 5, -1);
		
		
		$nCurrChr = 0;
		$nCurrLocusOnChr = 0;
		$nTotalLocus = 0;
		//echo("$nPopId, $nIndId\n");
		
		$nTotalCol = count($arrColumns);
		//echo($nTotalCol.PHP_EOL);
		for( $nCol = 0; $nCol < $nTotalCol; $nCol ++ )
		{
			//echo $nTempCounter.PHP_EOL;
			$nCell = $arrColumns[$nCol];
			
			if ($nCell == -1) {
				$nCurrChr++;
				$nCurrLocusOnChr = 1;
				$nMarkersToSampleCount =0;
				continue;
			}
			else {

					//echo($nLnCount);
					//die();
					//echo($arrMarkersToSample[$nCurrChr-1][$nMarkersToSampleCount]."\t".$nCurrLocusOnChr.PHP_EOL);
					if ($arrMarkersToSample[$nCurrChr-1][$nMarkersToSampleCount] != ($nCurrLocusOnChr - 1)) { //dont sample this locus
						$nCurrLocusOnChr++;
						//echo("exclude\n");
						continue;
					}
					
					
					//echo("sample\n");	
									
					if ($nLnCount == 0) {
						$arrAlleleCount_bigA_Pop1[$nTotalLocus] = 0;
						$arrAlleleCount_bigA_Pop2[$nTotalLocus] = 0;
						$arrAlleleCount_littleA_Pop1[$nTotalLocus] = 0;
						$arrAlleleCount_littleA_Pop2[$nTotalLocus]  = 0;
						$arrAlleleName[$nTotalLocus] = "chr".$nCurrChr."_".$nCurrLocusOnChr;
						//echo("line 0");
						//die();
					}
					
					
				$nMarkersToSampleCount++;
				
				$sAllele = substr($nCell , 0, 1);
				//echo $sAllele.PHP_EOL;
				if ($sAllele == "A") {

					
					if ($nPopId == 1) {
						//echo "la1".PHP_EOL;
						$arrAlleleCount_bigA_Pop1[$nTotalLocus]++;//isset($arrAlleleCount_bigA_Pop1[$nTotalLocus])?  ($arrAlleleCount_bigA_Pop1[$nTotalLocus] + 1) : 0;
					}
					if ($nPopId == 2) {
						//echo "la2".PHP_EOL;
						$arrAlleleCount_bigA_Pop2[$nTotalLocus]++;//= isset($arrAlleleCount_bigA_Pop2[$nTotalLocus])?  ($arrAlleleCount_bigA_Pop2[$nTotalLocus] + 1) : 0;
					}
				}
				
				if ($sAllele == "a") {
					if ($nPopId == 1) {
						//echo "la3".PHP_EOL;
						$arrAlleleCount_littleA_Pop1[$nTotalLocus]++ ;// = isset($arrAlleleCount_littleA_Pop1[$nTotalLocus])?  ($arrAlleleCount_littleA_Pop1[$nTotalLocus] + 1) : 0;
					}
					if ($nPopId == 2) {
						//echo "la4".PHP_EOL;
						$arrAlleleCount_littleA_Pop2[$nTotalLocus]++ ;// = isset($arrAlleleCount_littleA_Pop2[$nTotalLocus])?  ($arrAlleleCount_littleA_Pop2[$nTotalLocus] + 1) : 0;
					}
				}
				/*
				if ( $nLnCount == 0)
				{
				
					$arrAlleleName[$nTotalLocus] = "chr".$nCurrChr."_".$nCurrLocusOnChr;
					print_r($arrAlleleName);
					die();
				
				}
				else {
					print_r($arrAlleleName);
					die();
				}
				*/
				
				$nTotalLocus++;
				$nCurrLocusOnChr++;
				
				
			}
		}
	}
	
	else {
		$nCurrPop3++;
		
		if ($nCurrPop3 > $nPop3) {
			continue;
		}
		$arrColumns = array_slice(explode("\t", $sLn) , 5, -1);
		
		//echo($nTotalLocus.PHP_EOL);
		$nCurrChr = 0;
		$nCurrLocusOnChr = 0;
		$nTotalLocus = 0;
		$nTotalMarkers = count($arrAlleleName);
		
		if ($bFirstPop3) { // write genotypes.txt header
			 
			 
				fwrite( $fExGenotypesFile  , "id\t");
				//echo($nTotalMarkers);
				//die();
				for ($nMarker=0;$nMarker< $nTotalMarkers ;$nMarker++) {
					fwrite( $fExGenotypesFile  , $arrAlleleName[$nMarker]);
					if ($nMarker != $nTotalMarkers-1) {
						fwrite( $fExGenotypesFile  , "\t");
					}
				}
				//fwrite( $fExGenotypesFile  , PHP_EOL);
				$bFirstPop3=false;
		}
		
		//echo("$nPopId, $nIndId\n");
		$nTotalCol = count($arrColumns);
		if ($nChrSide == 1)
		{
			fwrite( $fExGenotypesFile  , PHP_EOL."$nIndId\t");
		}
		//echo($nTotalCol.PHP_EOL);
		for( $nCol = 0; $nCol < $nTotalCol; $nCol++ )
		{
			//echo $nTempCounter.PHP_EOL;
			$nCell = $arrColumns[$nCol];
			
			if ($nCell == "-1") {
				//echo("col$nCol\n" );
				$nCurrChr++;
				$nCurrLocusOnChr = 1;
				$nMarkersToSampleCount =0;
				continue;
			}
			else {
				/*
					if($nMarkersToSampleCount == 446 && $nCurrChr==1) {
						echo("chr1 sucess");
						die();
					}
					echo($nMarkersToSampleCount."\t".($arrMarkersToSample[$nCurrChr-1][$nMarkersToSampleCount])."\t".($nCurrLocusOnChr-1).PHP_EOL);
					if ($nCurrChr > 1) {
						die();
					}
			*/
					if ($arrMarkersToSample[$nCurrChr-1][$nMarkersToSampleCount] != ($nCurrLocusOnChr-1)) { //dont sample this locus
						$nCurrLocusOnChr++;
					
						continue;
					}
					
				$nMarkersToSampleCount++;
				
				$sAllele = substr($nCell , 0, 1);
				//echo $sAllele.PHP_EOL;
				if ($nChrSide == 0) {
					$arrAlleleTemp[$nTotalLocus] = ($sAllele=="A")? 1:2;//save this to temp array;
				}
				else {
				/*
					echo(count($arrAlleleTemp).PHP_EOL);
					die();
					*/
					fwrite( $fExGenotypesFile , "\"".$arrAlleleTemp[$nTotalLocus].",".(($sAllele=="A")? 1:2)."\"");
					
					//echo($nMarkersToSampleCount.PHP_EOL);
					
					if ($nMarkersToSampleCount != $nTotalMarkers-1) {
						fwrite( $fExGenotypesFile  , "\t");
						
					}
		
				}
				$nCurrLocusOnChr++;
				$nTotalLocus++;
			}
		}
		
		
		
	}
	
	$nChrSide = ($nChrSide==0)? 1:0;
	
 }
 
 //write priorallelefreq file:
 $nTotalMarkers = count($arrAlleleName);
 echo("Writing priorallelefreq...\n");
fwrite( $fPriorAlleleFreqFile , "LocusNames\tPop1\tPop2".PHP_EOL);
 for ($nMarker=0;$nMarker< $nTotalMarkers ;$nMarker++) {
	fwrite( $fPriorAlleleFreqFile , $arrAlleleName[$nMarker]."\t".($arrAlleleCount_bigA_Pop1[$nMarker]+0.5)."\t".($arrAlleleCount_bigA_Pop2[$nMarker]+0.5).PHP_EOL);
	fwrite( $fPriorAlleleFreqFile , $arrAlleleName[$nMarker]."\t".($arrAlleleCount_littleA_Pop1[$nMarker]+0.5)."\t".($arrAlleleCount_littleA_Pop2[$nMarker]+0.5).PHP_EOL);
 }

 fclose($fGenotypeFile);
 fclose($fPriorAlleleFreqFile);
fclose($fExGenotypesFile);

?>