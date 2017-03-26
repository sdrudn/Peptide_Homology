# Peptide_Homology

This tool can be used for searching homology between two sets of peptides or proteins.

# Dependencies:

Python 2.7

# What does this program do?
searches for homology between two sets of peptides or proteins.

# Protocol
A protein_library.txt file placed in your Python27 folder (C:\Python27). See example file.
A peptides.txt file placed in your Python27 folder (C:\Python27). See example file.

Put the Peptide_Homology.py file inside your Python 2.7 folder (e.g. C:\Python27)
Open command prompt on your computer and use the following commands.
Change directory to your Python27 folder by entering:
Cd C:\Python27
Conduct the homology search by entering ( homology can be between 1 and 100; the file names may be different in your case):

python.exe match_pep.py --proteins protein_library.txt --peptides bioactive_peptides.txt --pident_threshold 80 > 80_thres_similar_out.txt

This one uses amino acids with similar characters as the same. If you only want to use identical amino acids add --use_identical as parameter:

python.exe match_pep.py --proteins protein_library.txt --peptides bioactive_peptides.txt --use_identical --pident_threshold 80 > 80_thres_similar_out.txt
