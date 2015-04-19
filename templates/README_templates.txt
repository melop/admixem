1) open each template in excel, except for genes_template.txt, where an excel file is already provided to you.
2) fill in the desired information. 
IMPORTANT: For the genes file, please allow Excel to execute macros.
IMPORTANT: In the genes file, copy an existing row and paste to the end of the file, change values around and it will be automatically sorted for you.
3) Once you're happy with the file content, open a new text editor, copy and paste from excel into the text file.
3) If you are using a mac, check your files (using for example "less filename") for the presence of ^M symbol. If this is present, convert your files to the correct format using a program such as the mac2unix program:
mac2unix phenotypes_template.prn
IMPORTANT: sometimes excel will introduce quotation marks around strings of text it recognizes as formulas. You can see these using less ("less filename"). If these are present, they need to be removed, for example by:

perl -pi -e 's/"//g' filename
