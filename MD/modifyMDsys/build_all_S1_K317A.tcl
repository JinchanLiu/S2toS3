###################################### Warning ###################################
## This script requires vmd1.9.3 to run!!!!! 					##
## Later version of psfgen doesn't work for the complicated OEC patches!!!!	##
##################################################################################

# If you want to mutate a residue of PSII with 3wu2_PSII_S1_addwat1, you need to 
# separate that chain and write a new segment for that chain.
# Say, you want to mutate K317 of segname PROD to ALA, you need to build a new PROD,
# and combine it with everything else.
# However, you cannot separate things simply into seganme PROD and not segname PROD,
# because psfgen isn't smart enough - it doesn't apply all of the patches you have 
# had, and you have to reapply the patch to your system when you combine the part
# you modify (in this case segname PROA) with other parts of the system.
#
# With the way the patches are written to bond OEC and cofactor to the protein, we 
# really have to rebuild the whole thing to get it right.
#
# Please citethe following paper if you use this file:
#
#


package require psfgen
topology ../toppar/top_all36_prot.rtf
topology ../toppar/top_all36_lipid.rtf
topology ../toppar/top_all36_cgenff.rtf
topology ../toppar/top_all36_carb.rtf
topology ../toppar/toppar_ions_won.str
topology ../toppar/toppar_all36_prot_heme.str
topology ../toppar/toppar_iron_heme.str
topology ../toppar/top_all36_prot_tyrOX.rtf
topology ../toppar/toppar_PSII_OEC.str

# load 3wu2_PSII_S1_addwat1
mol new ../3wu2_PSII_S1_addwat1.psf 
mol addfile ../3wu2_PSII_S1_addwat1.pdb 

# write 3wu2_PSII_S1_addwat1 into different parts 
set protein [atomselect top "chain P and not segname PROA PROB PROC PROD"]
set PROA [atomselect top "segname PROA"]
set PROBC [atomselect top "segname PROB PROC"]
set PROD [atomselect top "segname PROD"]
set cofactor [atomselect top "chain C and not resid 54"]
set cofactor54 [atomselect top "chain C and resid 54"]
set lipidX [atomselect top "chain L"]
set IX [atomselect top "segname IX"]
set WX [atomselect top "segname WX WG"]
set others [atomselect top "not chain P C L and not segname IX WX WG OEC"]

mkdir buildS1MUT

$protein writepsf buildS1MUT/3wu2_PSII_S1_eq4_prot.psf; $protein writepdb buildS1MUT/3wu2_PSII_S1_eq4_prot.pdb
$PROA writepsf buildS1MUT/3wu2_PSII_S1_eq4_PROA.psf; $PROA writepdb buildS1MUT/3wu2_PSII_S1_eq4_PROA.pdb
$PROBC writepsf buildS1MUT/3wu2_PSII_S1_eq4_PROBC.psf; $PROBC writepdb buildS1MUT/3wu2_PSII_S1_eq4_PROBC.pdb
$PROD writepsf buildS1MUT/3wu2_PSII_S1_eq4_PROD.psf; $PROD writepdb buildS1MUT/3wu2_PSII_S1_eq4_PROD.pdb
$cofactor writepsf buildS1MUT/3wu2_PSII_S1_eq4_cofactor.psf; $cofactor writepdb buildS1MUT/3wu2_PSII_S1_eq4_cofactor.pdb
$cofactor54 writepsf buildS1MUT/3wu2_PSII_S1_eq4_cofactor54.psf; $cofactor54 writepdb buildS1MUT/3wu2_PSII_S1_eq4_cofactor54.pdb
$lipidX writepsf buildS1MUT/3wu2_PSII_S1_eq4_lipidX.psf; $lipidX writepdb buildS1MUT/3wu2_PSII_S1_eq4_lipidX.pdb
$IX writepsf buildS1MUT/3wu2_PSII_S1_eq4_IX.psf; $IX writepdb buildS1MUT/3wu2_PSII_S1_eq4_IX.pdb
$WX writepsf buildS1MUT/3wu2_PSII_S1_eq4_WX.psf; $WX writepdb buildS1MUT/3wu2_PSII_S1_eq4_WX.pdb
$others writepsf buildS1MUT/3wu2_PSII_S1_eq4_others.psf; $others writepdb buildS1MUT/3wu2_PSII_S1_eq4_others.pdb

# Now we start to build things
resetpsf

pdbalias residue HIS HSD
pdbalias residue HIE HSE
pdbalias residue HID HSD
pdbalias residue HIP HSP
pdbalias atom ILE CD1 CD

##################################################################################
# build the segname PROA and remember, we need to apply all of the patches again
segment PROA {
  pdb buildS1MUT/3wu2_PSII_S1_eq4_PROA.pdb
  first ACE
  last CTER
}
coordpdb buildS1MUT/3wu2_PSII_S1_eq4_PROA.pdb PROA
patch GLUP PROA:65
guesscoord
regenerate angles dihedrals

writepsf buildS1MUT/3wu2_PSII_S1_PROA.psf
writepdb buildS1MUT/3wu2_PSII_S1_PROA.guess.pdb

# psfgen is dumb that it doesn't read the position of the H added to E65 COO-, 
# rather it just guess the position of H, so we need to copy the coordinate of 
# that H from what we have already had.
mol new buildS1MUT/3wu2_PSII_S1_PROA.psf
mol addfile buildS1MUT/3wu2_PSII_S1_PROA.guess.pdb
set E65H [atomselect 1 "segname PROA and resid 65"]
mol new buildS1MUT/3wu2_PSII_S1_eq4_PROA.pdb
set ref [atomselect 2 "segname PROA and resid 65"]
puts [$E65H get {x y z}]
puts [$ref get {x y z}]
$E65H set {x y z} [$ref get {x y z}]
set all [atomselect 1 all]
$all writepdb buildS1MUT/3wu2_PSII_S1_PROA.pdb
# finish building segname PROA
##################################################################################

resetpsf
##################################################################################
# now we are going to build segname PROD
segment PROD {
  pdb buildS1MUT/3wu2_PSII_S1_eq4_PROD.pdb
  first NTER
  last CTER
  # Here is where we mutate things!!
  mutate 317 ALA
}
coordpdb buildS1MUT/3wu2_PSII_S1_eq4_PROD.pdb PROD
guesscoord
regenerate angles dihedrals
writepsf buildS1MUT/3wu2_PSII_S1_K317A_PROD.psf
writepdb buildS1MUT/3wu2_PSII_S1_K317A_PROD.pdb
# finish building segname PROD
##################################################################################



resetpsf
##################################################################################
# Now we are going to combine everything
readpsf S1_OEC_addangle.psf pdb S1_OEC_addangle.pdb
readpsf buildS1MUT/3wu2_PSII_S1_PROA.psf pdb buildS1MUT/3wu2_PSII_S1_PROA.pdb
readpsf buildS1MUT/3wu2_PSII_S1_eq4_PROBC.psf pdb buildS1MUT/3wu2_PSII_S1_eq4_PROBC.pdb
readpsf buildS1MUT/3wu2_PSII_S1_K317A_PROD.psf pdb buildS1MUT/3wu2_PSII_S1_K317A_PROD.pdb
readpsf buildS1MUT/3wu2_PSII_S1_eq4_prot.psf pdb buildS1MUT/3wu2_PSII_S1_eq4_prot.pdb
readpsf buildS1MUT/3wu2_PSII_S1_eq4_cofactor.psf pdb buildS1MUT/3wu2_PSII_S1_eq4_cofactor.pdb 
readpsf buildS1MUT/3wu2_PSII_S1_eq4_cofactor54.psf pdb buildS1MUT/3wu2_PSII_S1_eq4_cofactor54.pdb
readpsf buildS1MUT/3wu2_PSII_S1_eq4_lipidX.psf pdb buildS1MUT/3wu2_PSII_S1_eq4_lipidX.pdb
readpsf buildS1MUT/3wu2_PSII_S1_eq4_IX.psf pdb buildS1MUT/3wu2_PSII_S1_eq4_IX.pdb
readpsf buildS1MUT/3wu2_PSII_S1_eq4_WX.psf pdb buildS1MUT/3wu2_PSII_S1_eq4_WX.pdb
readpsf buildS1MUT/3wu2_PSII_S1_eq4_others.psf pdb buildS1MUT/3wu2_PSII_S1_eq4_others.pdb

# Don't forget to adjust charge!!!!!!!!!!!!!!!!!!!!!!
delatom ION 534

patch S1MN1 OEC:1 PROA:189 PROA:342 PROA:332
patch S1MN2 OEC:1 PROC:354 PROA:342 PROA:344
patch S1MN3 OEC:1 PROC:354 PROA:333
patch S1MN4 OEC:1 PROA:170 PROA:333
patch S1CA1 OEC:1 PROA:170 PROA:344 PROA:189


patch FEBI2 CX:53 PROA:272 PROA:215 PROD:268 PROD:214 CX:54
patch HEHI2 CX:51 PROE:23 PROF:24
patch HEHI2 CX:52 PROV:41 PROV:92


writepsf 3wu2_PSII_S1_K317A.psf
writepdb 3wu2_PSII_S1_K317A.pdb
# finish combining things together
##################################################################################


##################################################################################
# This is just checking that we have the correct charge
mol new 3wu2_PSII_S1_K317A.psf 
mol addfile 3wu2_PSII_S1_K317A.pdb

set all [atomselect top all]

puts [vecsum [$all get charge]]

set segnamelist [lsort -unique [$all get chain]]

foreach segname ${segnamelist} {
  set sel [atomselect top "chain $segname"]
  puts "$segname [vecsum [$sel get charge]]"
  $sel delete
}
# Finish checking
##################################################################################


quit



