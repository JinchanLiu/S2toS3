This folder in organized by Jinchan Liu (jinchan.liu09@gmail.com) on Oct 26th, 2023

It gives you an example of how you can use the 3wu2_PSII_S1_addwat1 system we built
to generate systems with mutations. In this case, we mutate the PROD-K317 to ALA.

Please cite the following paper if you use files in this folder:

 

The script only use relative path, and all files needed are included, so it should 
run wherever you download it. 

Use vmd1.9.3 to run it!!! The psfgen embedded in 1.9.4 is weird and it doesn't like 
the complicated patch I wrote for the OEC.

Here is the command to run it:

	vmd1.9.3 -dispdev text -e build_all_S1_K317A.tcl

If you do other mutations, DON'T FORGET TO ADJUST YOUR CHARGE!!!!!!!


