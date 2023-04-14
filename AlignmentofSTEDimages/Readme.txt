Readme file for rotation and alignment of STED images as performed in Jusyte et al., 2023

--- Properties of the aligned images ---
For the script to work the images that need to be aligned need to have the following properties:
- Images collected in two fluorescent channels
- size of the images 25x25 pixels, pixel size 20 nm (this is relevant for the plotting functions in the script)
- For side view images: For one of the channels a line is saved per AZs. All images of that AZ are aligned according to this line.
- For side view images: A direction line is saved per AZ. These direction lines indicate the location of the cytoplasm. 
- For top view images:  A circle line indicating the circular shape of the AZ.
- Images from all AZs combined in one stack and saved under the names:
	- for first and second channel: combStack-Direction_Genotype_Channel#_Larve(larve_range).tif (e.g. combStack-Side_Bdel_BRP_larve2-24.tif)
	- for alignment lines: Line_combStack-Direction_Genotype_larve(larve_range).tif (e.g. Line_combStack-Side_Bdel_BRP_Larve2-24.tif)
	- for direction lines: Direction_combStack_Direction_Genotype_Channel1_larve(larverange).tif (e.g. Direction_combStack-Side_Bdel_BRP_larve2-24.tif)
	- For circle lines: Circle_CombStack-Direction_Genotype_Channel#_larve(larverange).tif (e.g. Circle_CombStack-Top_CaM_BRP_Larve1-30.tif)
Different genotypes need to be saved in different folders with the respective genotype name.

--- Running and plot analysis of side view AZs --- 
Alignment of side view images and plotting can be performed with the script: RunandPlot_sideviewAnalysis.m
Make the following adjustments to the script:
- Line 4: Give the name of the folder path
- Line 8: enter name of first genotype
- Line 9-10: Enter name of the first and second fluorescent channel
- Line 13: Indicate the number of larve included in the stack for the first genotype as a list with the number of the first larve, and the last larve
- Line 18: enter the name of the second genotype
- Line 19: Same as in 13 but for the second genotype

--- Running and Plot analysis of top viewed AZs ---
Alignment of top view images and plotting can be performed with the script: RunandPlot_TopviewAnalysis.m
Make the following adjustments to the script:
- Line 3: Give the name of the folder path
- Line 8: enter name of first genotype
- Line 9-10: Enter name of the first and second fluorescent channel
- Line 12: Indicate the number of larve included in the stack for the first genotype as a list with the number of the first larve, and the last larve
- Line 26: enter the name of the second genotype
- Line 27: Same as in 12 but for the second genotype
- Line 128: enter name of first genotype
- Line 129: Indicate the number of larve included in the stack for the first genotype as a list with the number of the first larve, and the last larve
- Line 133: enter the name of the second genotype
- Line 134: Same as in 12 but for the second genotype

--- For bootstrapping of line profiles --- 
To bootstrap the line profiles (through the middle of the AZs) as described in Justyte et al., 2023 run the following code: Run_bootstraplines_topandside.m
This script runs both images obtained in the topview and in side view
Make the following adjustments to the script:
- Line 6: Give the name of the folder path
- Line 9: enter name of first genotype
- Line 10-11: Enter name of the first and second fluorescent channel
- Line 14: Indicate the number of larve included in the stack for the first genotype as a list with the number of the first larve, and the last larve
- Line 22: enter the name of the second genotype
- Line 23: Same as in 14 but for the second genotype
- Line 73: Give the name of the folder path
- Line 75: enter name of first genotype
- Line 76: Indicate the number of larve included in the stack for the first genotype as a list with the number of the first larve, and the last larve
- Line 80: enter the name of the second genotype
- Line 81: Same as in 76 but for the second genotype






