# Profex-.str-to-TOPAS-.str-converter
a python script for the batch conversion of profex .str to TOPAS format .str 

thanks co-pilot

# Purpose
Batch conversion of profex str files to a format readable by TOPAS, e.g. the structures provided when Profex is installed. Full structures from V5.5 were used for testing.
Every unique profex_format.str will produce a uniquely named TOPAS_format.str.
The converter rejects ambiguity in a profex str and will skip the problematic site, or emit a warning. Very good example: arsenolite. 
The As site is written in the profex str as Wyckoff e x 0.3529 y 0 z 0 for space group 227:1. This is ambiguous at best, wrong at worst. The converter will skip this site as other sources of information are required to resolve this ambiguity. 
Any human looking at the profex str would be rightly confused. 

# How to use
•	Ensure path for SPACEGRP.DAT from profex install is set correctly in the .py file

•	DEFAULT_SPACEGRP_DAT = r"C:\Program Files\Profex5\BGMNwin\SPACEGRP.DAT"

•	Open command prompt

•	Download and save .py

•	Download and install python 

•	In command prompt:

    cd C:\STRING1 where STRING is the directory for the .py file

    python ProfexStr_to_TOPAS_v2.py ^

    --in "C:\Program Files\Profex5\Structures\SUB-FOLDER” ^ where SUB-FOLDER is the sub-folder of profex str to convert in the main install.

    --out "C:\STRING2\Profex_Converted_Str_to_TOPAS\SUB-FOLDER where STRING2 is wherever you want the output to go and SUB-FOLDER is the same sub-folder in profex 
  

•	THE USER SHOULD VERIFY THE OUTPUT IS CORRECT AND MATCHES THE INTENT OF THE PROFEX STR. THERE MAY BE FRINGE CASES I HAVE MISSED. 


# Known failures to convert
•	Arsenolite.str

•	smectitedi2wfix1.str  turbostratic disorder / layer modelling. Out of scope, for now…

•	nontronite15a.str  turbostratic disorder / layer modelling. Out of scope, for now….

•	CRYOLITE.STR  Wyckoff a and b both given as 0, 0, 0. Impossible for 14 P12_1/c1

•	CSH-0625.str

•	CoFe2O4.str

•	Cu2O.str

•	ZnAl2O4.str









