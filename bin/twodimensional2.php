<?php
error_reporting(E_ALL & ~E_NOTICE);
ini_set('display_errors', '1');
ini_set('memory_limit', '9999M');

$sGenotypeFile  = "genotypes.txt";
$sPhenotypeFile = "outcomevars.txt";
$sOutFile = "twodimscan.csctl.txt";
$nPhenotypeFilter = 0;
$sPhenotypeCol = 0;
$nSampleSize = 200;
$nSampleEvery = 5;

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-s':
            $sGenotypeFile   = array_shift($argv);
            break;
        case '-o':
            $sOutFile   = array_shift($argv);
            break;
		case '-p':
			$sPhenotypeCol  = array_shift($argv);
			break;

    }
}
 
/* print them */
echo("This is php. $sGenotypeFile, $sOutFile, $sPhenotypeCol");

//read marker subset file:
$fPhenotypeFile = fopen($sPhenotypeFile, "r");
$arrIndividualsCasesControls = array();

$nLn = -1;
$nNumCases=0;
$nNumControls=0;
 
while(($sLn = fgets( $fPhenotypeFile  ))!==false) {
	if ($nLn==-1) {
		$nLn++;
		$arrIndividualsCasesControls[$nLn] = false; // don't sample the first line
		continue;
	}
	//echo($sLn);
	if (trim($sLn) == "" ) {
		continue;
	}
	$arrCols = explode("\t", $sLn);
	$nPhenoValue = $arrCols[$sPhenotypeCol];
	
	$arrIndividualsCasesControls[$nLn] = ($nPhenoValue == $nPhenotypeFilter)? true:false; // true then is case, false is control
	if ($nPhenoValue == $nPhenotypeFilter) {
		$nNumCases++;
	}
	else {
	    $nNumControls++;
	}
	
		$nLn++;
}

fclose($fPhenotypeFile );

//print_r($arrIndividualsCasesControls);
//die();
 $fGenotypeFile = fopen($sGenotypeFile , "r"); // source
 
 $nLnCount = -1;
 $sLn = "";
 $arrHeaders = array();
 $arrIndv = array();
 $nLoci = 0;
 
 
 
 while(($sLn = fgets( $fGenotypeFile  ))!==false) {
	 $nLnCount++;
	 
	 if ($nLnCount == 0) {//header line
		$arrHeaders = explode("\t", trim($sLn)); // save header
		array_shift($arrHeaders); // remove first column.
		//print_r($arrHeaders);
		//die();
		continue;
	 }
	 
	// echo(1);
	 //check if need to sample this individual is "case"
	 if (!isset($arrIndividualsCasesControls[$nLnCount])) {
		break;
	 }
	 
	// echo(2);
	 /*
	 if (!$arrIndividualsCasesControls[$nLnCount]) {
		continue;
	 }
	 */
	 
	// echo(3);
	if (trim($sLn) == "") {
		continue;
	}
	
	// echo("4 | ");

	$arrColumns = explode("\t", $sLn);
	
	//echo($sLn);
	$arrLoci = array();
	$arrColLen = count($arrColumns);
	$nLoci = $arrColLen - 1;
	
	//echo($arrColLen." ");
	
	
	for ($i=1; $i<$arrColLen; $i++) {
		$sVal = $arrColumns[$i];
		if (trim($sVal)=="") {
			$nLoci--;
			continue;
		}
		$arrVals = explode(",", str_replace("\"", "", $sVal));
		//print_r($arrVals[1]);
		
		$nVal = intval($arrVals[0]) + intval($arrVals[1]);
		//print_r($arrVals);
		//die();
		array_push($arrLoci , $nVal);
		//print_r($arrLoci);
	}
	
	array_push( $arrIndv , $arrLoci);
	//echo(count($arrIndv)." ");
 }
 

 //start two dimensional scan!
 //print_r($arrIndv);
 
 $arrRet = array();
 $nTotalInd= ($nSampleSize < count($arrIndv))? $nSampleSize: count($arrIndv);
 
  echo("\nIndividuals as cases: $nInd");
 

	 for ($i=0;$i< $nLoci ; $i+=$nSampleEvery) {
		$arrTemp1 = explode("_", $arrHeaders[$i]);
		$sChr = $arrTemp1[0];
	 
		for ($j=($i+1);$j < $nLoci ; $j+=$nSampleEvery) {
			$arrTemp2 = explode("_", $arrHeaders[$j]);
			$sThisChr = $arrTemp2[0];
			if ($sThisChr == $sChr ) {
				continue;
				}
			$sKey = $i."_".$j.":".$arrHeaders[$i]."~".$arrHeaders[$j];
			$arrRet[$sKey] = array();
	 }
}

 for ($nInd=0; $nInd < $nTotalInd; $nInd++) {
 
	echo("Ind: $nInd \n");
	 for ($i=0;$i< $nLoci - 1; $i+=$nSampleEvery) {
	 	$arrTemp1 = explode("_", $arrHeaders[$i]);
		$sChr = $arrTemp1[0];
		
		for ($j=($i+1);$j < $nLoci ; $j+=$nSampleEvery) {
			$arrTemp2 = explode("_", $arrHeaders[$j]);
			$sThisChr = $arrTemp2[0];
			if ($sThisChr == $sChr ) {
				continue;
				}
				
			$sKey = $i."_".$j.":".$arrHeaders[$i]."~".$arrHeaders[$j];
			
			$nDist = abs( $arrIndv[$nInd][$i] - $arrIndv[$nInd][$j]) ;//distance between loci, could be 2, 1 or 0
			//$nDist = ($arrIndv[$nInd][$i] == $arrIndv[$nInd][$j])? 1:0;
			if ( !isset($arrRet[$sKey][$nDist])) {
				$arrRet[$sKey][$nDist] = array();
				$arrRet[$sKey][$nDist][0] = 0;//num of cases 
				$arrRet[$sKey][$nDist][1] = 0;//num of controls
				//arrRet[$sKey][$nDist][2] = 0;//controls + cases
			}
			//$arrRet[$sKey][$nDist][2] += 1;
			if ($arrIndividualsCasesControls[$nInd]) { //case
				$arrRet[$sKey][$nDist][0] += 1;
			}
			else { //control
				$arrRet[$sKey][$nDist][1] += 1;
			}
		}
	 }
}

$fOutFile = fopen($sOutFile , "w");
	//fputs($fOutFile,  "PairName\t". "Dist=0" ."\t". "Dist=1" . "\t" . "Dist=2" . "\t". "ChiSq" . "\t". "zScore" . "\n");
	fputs($fOutFile,  "PairName\t". "ChiSq1" ."\t"."ChiSq2". "\tCountCase2\tCountCase1\tCountCase0\tCountControl2\tCountControl1\tCountControl0\t"."\n");
$nExpectedCase2 = 0.125;
$nExpectedCase1 = 0.5;
$nExpectedCase0 = 0.375;

$nExpectedControl2 = 0.125;
$nExpectedControl1 = 0.5;
$nExpectedControl0 = 0.375;

foreach ($arrRet as $sPairName => $arrDist) {

	//expected frequencies: 0.375:0, 0.5:1, 0.125:2
	//$nAvg = array_sum($arrDist)/3;
	/*
	$nCount0 = (isset($arrDist["0"])? $arrDist["0"]:0);
	$nCount1And2 = $nTotalInd - $nCount0;
	if ($nCount1And2 != 0) {
	$nCount1 = (isset($arrDist["1"])? $arrDist["1"]:0)/$nCount1And2;
	$nCount2 = (isset($arrDist["2"])? $arrDist["2"]:0)/$nCount1And2;
	//$nStd = sqrt((pow($nCount0 - $nAvg , 2) + pow($nCount1 - $nAvg , 2) + pow($nCount2 - $nAvg , 2))/3);
	//$nChiSq = pow($nCount0 - 0.375, 2)/0.375 + pow($nCount1 - 0.5, 2)/0.5 + pow($nCount2 - 0.125, 2)/0.125;
	$nChiSq = pow($nCount1 - 0.8, 2)/0.8 + pow($nCount2 - 0.2, 2)/0.2;
	}
	else {
		$nChiSq = 0;
	}
	*/
	/*
	$nCount1 = (isset($arrDist["0"])? $arrDist["0"]:0)/$nTotalInd;
	$nCount1 = (isset($arrDist["1"])? $arrDist["1"]:0)/$nTotalInd;
	$nCount2 = (isset($arrDist["2"])? $arrDist["2"]:0)/$nTotalInd;
	$nChiSq = pow($nCount0 - 0.375, 2)/0.375 + pow($nCount1 - 0.5, 2)/0.5 + pow($nCount2 - 0.125, 2)/0.125;
	$zScore = ($nCount1 + $nCount2 - 0.625)/0.625;
	*/
	
	$nCountCase2 = (isset($arrDist["2"]["0"])? $arrDist["2"]["0"]:0)/$nNumCases;
	$nCountCase1 = (isset($arrDist["1"]["0"])? $arrDist["1"]["0"]:0)/$nNumCases;
	$nCountCase0 = (isset($arrDist["0"]["0"])? $arrDist["0"]["0"]:0)/$nNumCases;
	$nCountControl2 = (isset($arrDist["2"]["1"])? $arrDist["2"]["1"]:0)/$nNumControls;
	$nCountControl1 = (isset($arrDist["1"]["1"])? $arrDist["1"]["1"]:0)/$nNumControls;
	$nCountControl0 = (isset($arrDist["0"]["1"])? $arrDist["0"]["1"]:0)/$nNumControls;
	
	$nChiSq = pow($nCountCase2 - $nExpectedCase2, 2)/$nExpectedCase2 + pow($nCountCase1 - $nExpectedCase1, 2)/$nExpectedCase1 + pow($nCountCase0 - $nExpectedCase0, 2)/$nExpectedCase0+ pow($nCountControl2 - $nExpectedControl2, 2)/$nExpectedControl2+ pow($nCountControl1 - $nExpectedControl1, 2)/$nExpectedControl1+ pow($nCountControl0 - $nExpectedControl0, 2)/$nExpectedControl0;
	$nChiSq2 =  pow($nCountCase2 - $nExpectedCase2, 2)/$nExpectedCase2 + pow($nCountCase1 - $nExpectedCase1, 2)/$nExpectedCase1 + pow($nCountCase0 - $nExpectedCase0, 2)/$nExpectedCase0;
	//fputs($fOutFile,  "$sPairName\t". $arrDist["0"] ."\t". $arrDist["1"]. "\t" . $arrDist["2"] . "\t". $nChiSq . "\t". $zScore . "\n");
	fputs($fOutFile,  "$sPairName\t". $nChiSq . "\t". $nChiSq2. "\t$nCountCase2\t$nCountCase1\t$nCountCase0\t$nCountControl2\t$nCountControl1\t$nCountControl0\t". "\n");
}

fclose($fOutFile);
?>