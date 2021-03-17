<?php
$sPhenotypeFile  = "";
$sOutcomeVarFile = "";

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-s':
            $sPhenotypeFile   = array_shift($argv);
            break;
        case '-o':
            $sOutcomeVarFile = array_shift($argv);
            break;

    }
}
 
/* print them */
echo("This is php. $sPhenotypeFile, $sOutcomeVarFile");
 

 $fPhenotypeFile = fopen($sPhenotypeFile , "r");
 $fOutcomeVarFile = fopen($sOutcomeVarFile , "w");
 $sLn = "";
 $bHeaderWritten = false;
 
 if ($fPhenotypeFile == false ) {
	echo("Input file not found\n");
 }
 

 
 while(($sLn = fgets( $fPhenotypeFile ))!==false) {
	$arrColumns = explode("\t", $sLn);
	//echo($sLn);
	
	if ($arrColumns[0] == "id") { // if header
	
		if (! $bHeaderWritten) { // write header
		
			for ($i=0;$i<5;$i++)
			{
				array_shift($arrColumns );
			}
			//print_r($arrColumns);
			//break;
			
			for ($i=0;$i<=count($arrColumns)-2;$i++) {
			
			
				fwrite($fOutcomeVarFile , "\"".$arrColumns[$i]."\"");
				if ($i != count($arrColumns)-2) {
					fwrite($fOutcomeVarFile , "\t");
				}
			}
			
			fwrite($fOutcomeVarFile, PHP_EOL);
			
			$bHeaderWritten = true;
			continue;
		}
		
		continue;
		
	}
	
	if ($arrColumns[1] != 3) { // if not pop 3, continue
		continue;
	}
	
	$arrPhenoVals= array_slice($arrColumns , 5, -1);
	fwrite($fOutcomeVarFile, implode("\t", $arrPhenoVals).PHP_EOL);
	
	
	
 }

 fclose($fPhenotypeFile);
 fclose( $fOutcomeVarFile);
 
 

?>
