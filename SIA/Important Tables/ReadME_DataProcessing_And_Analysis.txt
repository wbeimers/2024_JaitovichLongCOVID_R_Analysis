*I started R analysis by adding 3 new columns to output degreaser: 
1) MassError (TRUE if ppm Error is greater than 10)
2) FoundIn50Files (TRUE if unknown features are found in >= 167 samples (50% of lipid samples n=337),
			FALSE if unknown is found in less than 167 samples,
			NA if it is a known feature)
3)UniqueID column was added that includes "RT_mz" values

* From the Boxplot for abundance, i found out that the following samples had the lowest median abundance
and i had error notes during sample prep so i will exclude some of them. these samples:

X20240808_SIA_Batch2_36.raw
X20240729_SIA_Batch6_16.raw
X20240729_SIA_Batch6_58.raw
X20240801_SIA_Batch7_2.raw
X20240809_SIA_Batch1_1030.raw

*The following samples, I didn't have notes about them but had very low median abundance (i kept them):                   
X20240729_SIA_Batch6_51.raw                  
X20240806_SIA_Batch10_102.raw   
X20240809_SIA_Batch1_5.raw              
 

*Also, I re-injected some samples from Batch 1 and 2 and I will keep only the following ones 
for further analysis:

X20240809_SIA_Batch1_100.raw
X20240809_SIA_Batch1_1030.raw

_____________________________________

341 study samples extracted + 51 QC samples (7 raw files from study samples were excluded from analysis, 5 due to low ext. vol, 2 were replicate injections).
1730 quantifiable features (850 annotated, the rest is unidentified)
 
