cd $ENIGMADIR/shellscripts/NumberOfRelatedSpecies			#For qsub, so we get back to this directory
export PATH=$ENIGMADIR/bin:$PATH							#So we can find the executables

rm accuracy_1
rm accuracy_2
rm accuracy_3
rm accuracy_4
rm accuracy_5
rm accuracy_6
rm accuracy_7
rm accuracy_8

for ((i=0;i<1;i++))
do

#generate 2 node tree evidence
rm temp																			#Generator won't overwrite data without permission
GeneratorL1_constgenes inputL1_25genes.gene temp > genout_L1.stk				#Generate L1UTR evidence
L1_2_stadn genout_L1.stk > genout_L1_trans.stk									#Add transitions to L1 evidence
strandrandomly_L1 genout_L1_trans.stk > genout_L1_trans_strands.stk				#Switch some genes to reverse strand
L1stk2L0stk genout_L1_trans_strands.stk > genout_newL0_trans_strands.stk
newL0stk2oldL0stk genout_newL0_trans_strands.stk > genout_L0_trans_strands.stk
getevidence genout_L0_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1,((((E1:.1)E,(F1:.1)F)EF:.1,((G1:.1)G:.1,(H1:.1)H:.1)GH:.1)EFGH:.1)ABCDEFGH:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 7 B C D E F G H .6 .6 .6 .6 .6 .6 .6 > species8.stk
getevidence genout_L0_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1,((((E1:.1)E,(F1:.1)F)EF:.1,((G1:.1)G:.1)GH:.1)EFGH:.1)ABCDEFGH:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 6 B C D E F G .6 .6 .6 .6 .6 .6 .6 > species7.stk
getevidence genout_L0_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1,((((E1:.1)E,(F1:.1)F)EF:.1)EFGH:.1)ABCDEFGH:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 5 B C D E F .6 .6 .6 .6 .6 > species6.stk
getevidence genout_L0_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1,((((E1:.1)E)EF:.1)EFGH:.1)ABCDEFGH:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 4 B C D E .6 .6 .6 .6 > species5.stk
getevidence genout_L0_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 3 B C D .6 .6 .6 > species4.stk
getevidence genout_L0_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1)CD:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 2 B C .6 .6 > species3.stk
getevidence genout_L0_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 1 B .6 > species2.stk
getevidence genout_L0_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 0 > species1.stk

xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg species1.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg species1.stk -ar > arout													#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L0_trans_strands.stk >> accuracy_1							#report accuracy


xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg species2.stk -x expanded.eg		#expand grammar so we can transpose some branches
maketranspose expanded.eg species backwards AB > transpose.eg									#transpose the branch ending AB
xrate -g transpose.eg species2.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species2.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L0_trans_strands.stk >> accuracy_2							#report accuracy

xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg species3.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species backwards AB ABCD > transpose.eg							#transpose the branch ending AB
xrate -g transpose.eg species3.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species3.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L0_trans_strands.stk >> accuracy_3							#report accuracy

xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg species4.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species backwards AB ABCD > transpose.eg							#transpose the branch ending AB
xrate -g transpose.eg species4.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species4.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L0_trans_strands.stk >> accuracy_4							#report accuracy

xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg species5.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species backwards AB ABCD ABCDEFGH > transpose.eg					#transpose the branch ending AB
xrate -g transpose.eg species5.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species5.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L0_trans_strands.stk >> accuracy_5							#report accuracy

xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg species6.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species backwards AB ABCD ABCDEFGH > transpose.eg					#transpose the branch ending AB
xrate -g transpose.eg species6.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species6.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L0_trans_strands.stk >> accuracy_6							#report accuracy

xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg species7.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species backwards AB ABCD ABCDEFGH > transpose.eg					#transpose the branch ending AB
xrate -g transpose.eg species7.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species7.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L0_trans_strands.stk >> accuracy_7							#report accuracy

xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg species8.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species backwards AB ABCD ABCDEFGH > transpose.eg					#transpose the branch ending AB
xrate -g transpose.eg species8.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species8.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L0_trans_strands.stk >> accuracy_8							#report accuracy

done

