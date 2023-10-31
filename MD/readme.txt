This folder is organized by Jinchan Liu (jinchan.liu09@gmail.com) on Oct 26th, 2023

It includes the files we used to run the Molecular Dynamics (MD) of the replenishment
of W10, after W10 binds to Mn4 during the S2 to S3 transiton.

Please cite the following paper if you use files in this folder:


 

##################################################################################
## prarameter and topology files						##
##################################################################################
toppar/*

All prarameter and topology files we used can be found in the toppar folder.
Note: Because the atom names in the BCR and lipids solved by 3wu2 don't match with
charmm36, we wrote the topologies for them to include these cofactors in our system. 
Anyone who needs a PSII system should include them too - if they are solved in an Xray
structure, they are stable enough to serve some function!
##################################################################################


##################################################################################
## structure and coordinate files						##
##################################################################################
3wu2_PSII_S1_addwat1.p*
3wu2_PSII_S2DW.p*
3wu2_PSII_S2DW_EQed.pdb

We included the system we build at S1nYz (S1) and 
				   S2+cYzdot with W10 binding to Mn4 (S2DW)

The system we build includes the following parts:
chain O: the OEC (segname OEC)
chain P: the amino acids (segname PROA to PROZ)
chain C: cofactors solved by 3wu2 (segname CX) and 4v62 (segname CY)
chain L: lipids solved by 3wu2 (segname LX) and 4v62 (segname LY)
chain I: ions solved by 3wu2 (segname IX) and neutralizing the system (segname ION)
chain W: water solved by 3wu2 (segname WX), added by us (segname WA) and solvating the system (segname WT*)
chain G: glycolipids (segname G**) by charmmgui
chain M: other lipids by charmmgui
See BuildSys.docx and the Supplement of our paper for more details.
##################################################################################


##################################################################################
## namd files									##
##################################################################################
3wu2_PSII_S2DW_SW_rep1.conf
S1_restrain.in
S2DW_restrain.in
3wu2_PSII_S2DW_eq5.restart.xsc

An example of the namd configuration file we used to sample water is included 
*.in are the special restraints we applied.
3wu2_PSII_S2DW_eq5.restart.xsc includes the periodic boundary from previous eq step.
(Although in the Supplement I wrote 6 steps - the 6th step is eq5 - I didn't make a mistake here)
##################################################################################



