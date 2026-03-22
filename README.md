# Profex.str to TOPAS.str Converter
a python script for the batch conversion of profex .str to TOPAS format .str 

thanks co-pilot

[![DOI](https://zenodo.org/badge/1187204487.svg)](https://doi.org/10.5281/zenodo.19134548)

# Purpose
Batch conversion of profex str files to a format readable by TOPAS, e.g. the structures provided when Profex is installed. Full structures from V5.5 were used for testing.
Every unique profex_format.str will produce a uniquely named TOPAS_format.str.
The converter rejects ambiguity in a profex str and will skip the problematic site, or emit a warning. Very good example: arsenolite. 
The As site is written in the profex str as Wyckoff e x 0.3529 y 0 z 0 for space group 227:1. This is ambiguous at best, wrong at worst. The converter will skip this site as other sources of information are required to resolve this ambiguity. 
Any human looking at the profex str would be rightly confused. 

# How to use

*Ensure Profex is installed*

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
  

# THE USER SHOULD VERIFY THE OUTPUT IS CORRECT AND MATCHES THE INTENT OF THE PROFEX STR. THERE MAY BE FRINGE CASES I HAVE MISSED. 
Most issues during testing were the result of formatting differences in the profex str or typos. For example, SpacegropuNo instead of SpacegroupNo, or spaces after a value where EVERY other .str has no space. 
The original profex data is preserved as a comment for ALL atomic sites for this reason. 


# Known failures to convert
•	Arsenolite.str

•	smectitedi2wfix1.str    which has turbostratic disorder / layer modelling. Out of scope, for now…

•	nontronite15a.str     which has turbostratic disorder / layer modelling. Out of scope, for now….

•	CRYOLITE.STR  Wyckoff a and b both given as 0, 0, 0 in the .str file. This is impossible for 14 P12_1/c1 

•	CSH-0625.str

•	CoFe2O4.str

•	Cu2O.str

•	ZnAl2O4.str

# Conversion

As much of the original profex.str is retained as a comment for the user to confirm (e.g. numeric space group number:setting or HM notation). Phase names generated uniquely so no 2 files within the same batch conversion will be the same. Parameters derived from phase names.

Coordinates output as mod1 to avoid negatives

Beq = 100 * TDS

Beq 1 inserted when no TDS given as this is the TOPAS gui default for Beq handling when a structure or cif is imported where none are defined. Note, if you copy such str into a text editor for launch mode, TOPAS will use Beq 0 if you don't explicity add a Beq value.

Profex min and max limits on unit cell applied. If none given +/- 1% default applied 

Fixed unit cell parameters by crystal system do not have generated parameter names or are constrained (e.g. a = b = c for cubic with same parameter name)

MVW(), CS_L and scale inserted with derived parameter names

Occupancy parameters are inserted as well when necessary and also have unique names to avoid constraint. 

Coordinates converted to fractions. SPACEGRP.DAT used to insert missing coordinates for special positions. Free coordinates are never guessed - the site is skipped if this is the case in the original .str












