<?php
$arrOpts = getopt("d:o:p:s:c:l:u:t:"); // d- dir, o - outfile, p - population, s- sex filter 0/1, c - column to sum, l - range lower bound, u - range upper bound, t - range step
$sFiles = $arrOpts["d"]."/*_phenotypes.txt";
$sOutFile = $arrOpts["o"];
$nPopFilter = $arrOpts["p"];
$nSexFilter = $arrOpts["s"];
$nColToSum = $arrOpts["c"];
$arrValues = range($arrOpts["l"], $arrOpts["u"], $arrOpts["t"]); // phenotype values 0 - 16
$arrAllData = array();
$hOutFile = fopen($sOutFile , "w");

foreach (glob($sFiles) as $sPhenoFile) {
    //echo($sPhenoFile.PHP_EOL);
    
    preg_match_all('!\d+!', basename($sPhenoFile), $matches);
    //print_r($matches);
	$nGen = $matches[0][0];
	$arrCounts = array();
	foreach($arrValues as $nVal) {
		$arrCounts[$nVal] = 0; //array_fill_keys($arrValues , 0);
	}
	$hPhenoFile = fopen( $sPhenoFile, "r");
	
	while (($sLn = fgets($hPhenoFile)) !== false) {
	
		if (trim($sLn) == "") {
			continue;
		}
		
		$arrFields = explode("\t" , trim($sLn));
		if ( $arrFields == "id" ) {
			continue;
		}
		
		if ( $arrFields[1] != $nPopFilter ) {
			continue;
		}
		
		if ( $arrFields[2] !=  $nSexFilter ) {
			continue;
		}
		
		$arrCounts[$arrFields[$nColToSum]] += 1;
	}
	
	$arrAllData[ $nGen  ] = $arrCounts;
	
	fclose($hPhenoFile);
	
	
}

ksort($arrAllData);
fwrite($hOutFile , "Gen");
$bHeaderWritten = false;

foreach($arrAllData as $nGen => $arrCounts) {
	if (!$bHeaderWritten) {
		foreach( $arrCounts as $nVal => $nCount) {
			fwrite($hOutFile , "\t".$nVal);
		}
		fwrite($hOutFile , PHP_EOL);
		$bHeaderWritten = true;
	}
	
	fwrite($hOutFile ,  $nGen);
	
	foreach( $arrCounts as $nVal => $nCount) {
		fwrite($hOutFile , "\t".$nCount);
	}
	
	fwrite($hOutFile , PHP_EOL);
}

fclose($hOutFile);
?>