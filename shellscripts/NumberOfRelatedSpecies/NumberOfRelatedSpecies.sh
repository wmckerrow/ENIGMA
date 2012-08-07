#As yet untested. I gotten eval installed on fir yet.

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
rm exonsensitivity_1
rm exonsensitivity_2
rm exonsensitivity_3
rm exonsensitivity_4
rm exonsensitivity_5
rm exonsensitivity_6
rm exonsensitivity_7
rm exonsensitivity_8

for ((i=0;i<10;i++))
do

#generate 2 node tree evidence
rm temp																			#Generator won't overwrite data without permission
Generator_L1UTR_constgenes input8.gene temp > genout_L1utr.stk						#Generate L1UTR evidence
L1utr_2_stadn genout_L1utr.stk > genout_L1utr_trans.stk								#Add transitions to L1utr evidence
strandrandomly_L1utr genout_L1utr_trans.stk > genout_L1utr_trans_strands.stk			#Switch some genes to reverse strand
getevidence genout_L1utr_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1,((((E1:.1)E,(F1:.1)F)EF:.1,((G1:.1)G:.1,(H1:.1)H:.1)GH:.1)EFGH:.1)ABCDEFGH:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 7 B C D E F G H .5 .5 .5 .5 .5 .5 .5 > species8.stk
getevidence genout_L1utr_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1,((((E1:.1)E,(F1:.1)F)EF:.1,((G1:.1)G:.1)GH:.1)EFGH:.1)ABCDEFGH:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 6 B C D E F G .5 .5 .5 .5 .5 .5 .5 > species7.stk
getevidence genout_L1utr_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1,((((E1:.1)E,(F1:.1)F)EF:.1)EFGH:.1)ABCDEFGH:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 5 B C D E F .5 .5 .5 .5 .5 > species6.stk
getevidence genout_L1utr_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1,((((E1:.1)E)EF:.1)EFGH:.1)ABCDEFGH:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 4 B C D E .5 .5 .5 .5 > species5.stk
getevidence genout_L1utr_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1,(D1:.1)D:.1)CD:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 3 B C D .5 .5 .5 > species4.stk
getevidence genout_L1utr_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1,(((C1:.1)C:.1)CD:.1)ABCD:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 2 B C .5 .5 > species3.stk
getevidence genout_L1utr_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1,((B1:.1)B:.1)AB:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 1 B .5 > species2.stk
getevidence genout_L1utr_trans_strands.stk "#=GF NH (A1:.1,A2:.1,A3:.1,A4:.1)A;" 4 A1 A2 A3 A4 1 1 1 1 0 > species1.stk
GeneratePositions genout_L1utr_trans_strands.stk > pos
stk2gff_L1utr 0 genout_L1utr_trans_strands.stk genout pos
gff3_to_gtf.pl genout_A.gff > genout_A.gtf

xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg species1.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg species1backwards -ar > arout													#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L1utr_trans_strands.stk >> accuracy_1							#report accuracy
stk2gff_L1utr 0 xrateout.stk xrateout_1 pos													#convert stk to gff then to gtf for eval
gff3_to_gtf.pl xrateout_1_A.gff > xrateout_1_A.gtf
evaluate_gtf.pl genout_A.gtf xrateout_1_A.gtf > evalout_1									#run eval
parseeval evalout_1 exonsensitivity_1 79 2 79 3 11 3 12 3 > garbage							#extract exonsensitivity from eval output

xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg species2.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species reverse AB > transpose.eg									#transpose the branch ending AB
xrate -g transpose.eg species2.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species2.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L1utr_trans_strands.stk >> accuracy_2							#report accuracy
stk2gff_L1utr 0 xrateout.stk xrateout_2 pos													#convert stk to gff then to gtf for eval
gff3_to_gtf.pl xrateout_2_A.gff > xrateout_2_A.gtf
evaluate_gtf.pl genout_A.gtf xrateout_2_A.gtf > evalout_2									#run eval
parseeval evalout_2 exonsensitivity_2 79 2 79 3 11 3 12 3 > garbage							#extract exonsensitivity from eval output

xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg species3.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species reverse AB ABCD > transpose.eg							#transpose the branch ending AB
xrate -g transpose.eg species3.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species3.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L1utr_trans_strands.stk >> accuracy_3							#report accuracy
stk2gff_L1utr 0 xrateout.stk xrateout_3 pos													#convert stk to gff then to gtf for eval
gff3_to_gtf.pl xrateout_3_A.gff > xrateout_3_A.gtf
evaluate_gtf.pl genout_A.gtf xrateout_3_A.gtf > evalout_3									#run eval
parseeval evalout_3 exonsensitivity_3 79 2 79 3 11 3 12 3 > garbage							#extract exonsensitivity from eval output

xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg species4.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species reverse AB ABCD > transpose.eg							#transpose the branch ending AB
xrate -g transpose.eg species4.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species4.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L1utr_trans_strands.stk >> accuracy_4							#report accuracy
stk2gff_L1utr 0 xrateout.stk xrateout_4 pos													#convert stk to gff then to gtf for eval
gff3_to_gtf.pl xrateout_4_A.gff > xrateout_4_A.gtf
evaluate_gtf.pl genout_A.gtf xrateout_4_A.gtf > evalout_4									#run eval
parseeval evalout_4 exonsensitivity_4 79 2 79 3 11 3 12 3 > garbage							#extract exonsensitivity from eval output

xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg species5.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species reverse AB ABCD ABCDEFGH > transpose.eg					#transpose the branch ending AB
xrate -g transpose.eg species5.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species5.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L1utr_trans_strands.stk >> accuracy_5							#report accuracy
stk2gff_L1utr 0 xrateout.stk xrateout_5 pos													#convert stk to gff then to gtf for eval
gff3_to_gtf.pl xrateout_5_A.gff > xrateout_5_A.gtf
evaluate_gtf.pl genout_A.gtf xrateout_5_A.gtf > evalout_5									#run eval
parseeval evalout_5 exonsensitivity_5 79 2 79 3 11 3 12 3 > garbage							#extract exonsensitivity from eval output

xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg species6.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species reverse AB ABCD ABCDEFGH > transpose.eg					#transpose the branch ending AB
xrate -g transpose.eg species6.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species6.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L1utr_trans_strands.stk >> accuracy_6							#report accuracy
stk2gff_L1utr 0 xrateout.stk xrateout_6 pos													#convert stk to gff then to gtf for eval
gff3_to_gtf.pl xrateout_6_A.gff > xrateout_6_A.gtf
evaluate_gtf.pl genout_A.gtf xrateout_6_A.gtf > evalout_6									#run eval
parseeval evalout_6 exonsensitivity_6 79 2 79 3 11 3 12 3 > garbage							#extract exonsensitivity from eval output

xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg species7.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species reverse AB ABCD ABCDEFGH > transpose.eg					#transpose the branch ending AB
xrate -g transpose.eg species7.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species7.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L1utr_trans_strands.stk >> accuracy_7							#report accuracy
stk2gff_L1utr 0 xrateout.stk xrateout_7 pos													#convert stk to gff then to gtf for eval
gff3_to_gtf.pl xrateout_7_A.gff > xrateout_7_A.gtf
evaluate_gtf.pl genout_A.gtf xrateout_7_A.gtf > evalout_7									#run eval
parseeval evalout_7 exonsensitivity_7 79 2 79 3 11 3 12 3 > garbage							#extract exonsensitivity from eval output

xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg species8.stk -x expanded.eg  #expand grammar so we can transpose some branches
maketranspose expanded.eg species reverse AB ABCD ABCDEFGH > transpose.eg					#transpose the branch ending AB
xrate -g transpose.eg species8.stk -t trained.eg											#train the grammar parameters
xrate -g trained.eg species8.stk -ar > arout												#ancestral reconstruction
simplify_ar arout > xrateout																#create an stk with just the ar info
accuracy_L1 A xrateout genout_L1utr_trans_strands.stk >> accuracy_8							#report accuracy
stk2gff_L1utr 0 xrateout.stk xrateout_8 pos													#convert stk to gff then to gtf for eval
gff3_to_gtf.pl xrateout_8_A.gff > xrateout_8_A.gtf
evaluate_gtf.pl genout_A.gtf xrateout_8_A.gtf > evalout_8									#run eval
parseeval evalout_8 exonsensitivity_8 79 2 79 3 11 3 12 3 > garbage							#extract exonsensitivity from eval output

done

