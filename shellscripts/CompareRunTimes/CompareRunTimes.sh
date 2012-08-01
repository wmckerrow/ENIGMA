cd $ENIGMADIR/shellscripts/CompareRunTimes					#For qsub, so we get back to this directory
export PATH=$ENIGMADIR/bin:$PATH							#So we can find the executables

for ((i=0;i<5;i++))
do

rm times_L0
rm times_L0_trans
rm times_L0_trans_strands
rm times_L1_trans_strands
rm times_L1utr_trans_strands

#get 25 gene evidence
rm evidence_L1																	#Generator won't overwrite data without permission
GeneratorL1_constgenes inputL1_25genes.gene evidence_L1-2.stk > genout_L1.stk	#Generate L1 evidence
retreestk evidence_L1-2.stk "#=GF NH (((M1:.1,M2:.1,M3:.1,M4:.1)M:.1,(S1:.1)S:.1)MS:.1,(A1:.1)A:.1)AMS;" > evidence_L1.stk #We need dummy leaves so that xrate doesn't treat A and S as evidence
L1stk2L0stk evidence_L1.stk > evidence_L0.stk									#Convert to L0 evidence
eix2stadn evidence_L0.stk > evidence_L0_trans.stk								#Add transitions to L0 evidence
L1_2_stadn evidence_L1.stk > evidence_L1_trans.stk								#Add transitions to L1 evidence
strandrandomly_L1 evidence_L1_trans.stk > evidence_L1_trans_strands.stk			#Switch some genes to reverse strand
L1stk2L0stk evidence_L1_trans_strands.stk > evidence_newL0_trans_strands.stk	#So we have the same stranding
newL0stk2oldL0stk evidence_newL0_trans_strands.stk > evidence_L0_trans_strands.stk
L1stk2L1UTRstk evidence_L1_trans_strands.stk > evidence_L1utr_trans_strands.stk	#y becomes X

date >> times_L0																#start time
xrate -g $ENIGMADIR/grammars/L0.eg evidence_L0.stk -t trained.eg				#train the grammar parameters
xrate -g trained.eg evidence_L0.stk -ar > garbage								#ancestral reconstruction
date >> times_L0																#end time

date >> times_L0_trans															#start time
xrate -g $ENIGMADIR/grammars/L0_trans.eg evidence_L0_trans.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0_trans.stk -ar > garbage							#ancestral reconstruction
date >> times_L0_trans															#end time

date >> times_L0_trans_strands													#start time
xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg evidence_L0_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0_trans_strands.stk -ar > garbage					#ancestral reconstruction
date >> times_L0_trans_strands													#end time

date >> times_L1_trans_strands													#start time
xrate -g $ENIGMADIR/grammars/L1_trans_strands_4transpose.eg evidence_L1_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L1_trans_strands.stk -ar > garbage					#ancestral reconstruction
date >> times_L1_trans_strands													#end time

date >> times_L1utr_trans_strands												#start time
xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg evidence_L1utr_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L1utr_trans_strands.stk -ar > garbage				#ancestral reconstruction
date >> times_L1utr_trans_strands												#end time




#get 50 gene evidence
rm evidence_L1																	#Generator won't overwrite data without permission
GeneratorL1_constgenes inputL1_50genes.gene evidence_L1-2.stk > genout_L1.stk		#Generate L1 evidence
retreestk evidence_L1-2.stk "#=GF NH (((M1:.1,M2:.1,M3:.1,M4:.1)M:.1,(S1:.1)S:.1)MS:.1,(A1:.1)A:.1)AMS;" > evidence_L1.stk #We need dummy leaves so that xrate doesn't treat A and S as evidence
L1stk2L0stk evidence_L1.stk > evidence_L0.stk									#Convert to L0 evidence
eix2stadn evidence_L0.stk > evidence_L0_trans.stk								#Add transitions to L0 evidence
L1_2_stadn evidence_L1.stk > evidence_L1_trans.stk								#Add transitions to L1 evidence
strandrandomly_L1 evidence_L1_trans.stk > evidence_L1_trans_strands.stk			#Switch some genes to reverse strand
L1stk2L0stk evidence_L1_trans_strands.stk > evidence_newL0_trans_strands.stk	#So we have the same stranding
newL0stk2oldL0stk evidence_newL0_trans_strands.stk > evidence_L0_trans_strands.stk
L1stk2L1UTRstk evidence_L1_trans_strands.stk > evidence_L1utr_trans_strands.stk	#y becomes X

date >> times_L0													#start time
xrate -g $ENIGMADIR/grammars/L0.eg evidence_L0.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0.stk -ar > garbage					#ancestral reconstruction
date >> times_L0													#end time

date >> times_L0_trans															#start time
xrate -g $ENIGMADIR/grammars/L0_trans.eg evidence_L0_trans.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0_trans.stk -ar > garbage							#ancestral reconstruction
date >> times_L0_trans															#end time

date >> times_L0_trans_strands																			#start time
xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg evidence_L0_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0_trans_strands.stk -ar > garbage											#ancestral reconstruction
date >> times_L0_trans_strands																			#end time

date >> times_L1_trans_strands																			#start time
xrate -g $ENIGMADIR/grammars/L1_trans_strands_4transpose.eg evidence_L1_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L1_trans_strands.stk -ar > garbage											#ancestral reconstruction
date >> times_L1_trans_strands																			#end time

date >> times_L1utr_trans_strands																				#start time
xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg evidence_L1utr_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L1utr_trans_strands.stk -ar > garbage												#ancestral reconstruction
date >> times_L1utr_trans_strands																				#end time



#get 75 gene evidence
rm evidence_L1																	#Generator won't overwrite data without permission
GeneratorL1_constgenes inputL1_75genes.gene evidence_L1-2.stk > genout_L1.stk		#Generate L1 evidence
retreestk evidence_L1-2.stk "#=GF NH (((M1:.1,M2:.1,M3:.1,M4:.1)M:.1,(S1:.1)S:.1)MS:.1,(A1:.1)A:.1)AMS;" > evidence_L1.stk #We need dummy leaves so that xrate doesn't treat A and S as evidence
L1stk2L0stk evidence_L1.stk > evidence_L0.stk									#Convert to L0 evidence
eix2stadn evidence_L0.stk > evidence_L0_trans.stk								#Add transitions to L0 evidence
L1_2_stadn evidence_L1.stk > evidence_L1_trans.stk								#Add transitions to L1 evidence
strandrandomly_L1 evidence_L1_trans.stk > evidence_L1_trans_strands.stk			#Switch some genes to reverse strand
L1stk2L0stk evidence_L1_trans_strands.stk > evidence_newL0_trans_strands.stk	#So we have the same stranding
newL0stk2oldL0stk evidence_newL0_trans_strands.stk > evidence_L0_trans_strands.stk
L1stk2L1UTRstk evidence_L1_trans_strands.stk > evidence_L1utr_trans_strands.stk	#y becomes X

date >> times_L0													#start time
xrate -g $ENIGMADIR/grammars/L0.eg evidence_L0.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0.stk -ar > garbage					#ancestral reconstruction
date >> times_L0													#end time

date >> times_L0_trans															#start time
xrate -g $ENIGMADIR/grammars/L0_trans.eg evidence_L0_trans.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0_trans.stk -ar > garbage							#ancestral reconstruction
date >> times_L0_trans															#end time

date >> times_L0_trans_strands																			#start time
xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg evidence_L0_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0_trans_strands.stk -ar > garbage											#ancestral reconstruction
date >> times_L0_trans_strands																			#end time

date >> times_L1_trans_strands																			#start time
xrate -g $ENIGMADIR/grammars/L1_trans_strands_4transpose.eg evidence_L1_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L1_trans_strands.stk -ar > garbage											#ancestral reconstruction
date >> times_L1_trans_strands																			#end time

date >> times_L1utr_trans_strands																				#start time
xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg evidence_L1utr_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L1utr_trans_strands.stk -ar > garbage												#ancestral reconstruction
date >> times_L1utr_trans_strands																				#end time



#get 100 gene evidence
rm evidence_L1																	#Generator won't overwrite data without permission
GeneratorL1_constgenes inputL1_100genes.gene evidence_L1-2.stk > genout_L1.stk	#Generate L1 evidence
retreestk evidence_L1-2.stk "#=GF NH (((M1:.1,M2:.1,M3:.1,M4:.1)M:.1,(S1:.1)S:.1)MS:.1,(A1:.1)A:.1)AMS;" > evidence_L1.stk #We need dummy leaves so that xrate doesn't treat A and S as evidence
L1stk2L0stk evidence_L1.stk > evidence_L0.stk									#Convert to L0 evidence
eix2stadn evidence_L0.stk > evidence_L0_trans.stk								#Add transitions to L0 evidence
L1_2_stadn evidence_L1.stk > evidence_L1_trans.stk								#Add transitions to L1 evidence
strandrandomly_L1 evidence_L1_trans.stk > evidence_L1_trans_strands.stk			#Switch some genes to reverse strand
L1stk2L0stk evidence_L1_trans_strands.stk > evidence_newL0_trans_strands.stk	#So we have the same stranding
newL0stk2oldL0stk evidence_newL0_trans_strands.stk > evidence_L0_trans_strands.stk
L1stk2L1UTRstk evidence_L1_trans_strands.stk > evidence_L1utr_trans_strands.stk	#y becomes X

date >> times_L0													#start time
xrate -g $ENIGMADIR/grammars/L0.eg evidence_L0.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0.stk -ar > garbage					#ancestral reconstruction
date >> times_L0													#end time

date >> times_L0_trans															#start time
xrate -g $ENIGMADIR/grammars/L0_trans.eg evidence_L0_trans.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0_trans.stk -ar > garbage							#ancestral reconstruction
date >> times_L0_trans															#end time

date >> times_L0_trans_strands																			#start time
xrate -g $ENIGMADIR/grammars/L0_trans_strands_4transpose.eg evidence_L0_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L0_trans_strands.stk -ar > garbage											#ancestral reconstruction
date >> times_L0_trans_strands																			#end time

date >> times_L1_trans_strands																			#start time
xrate -g $ENIGMADIR/grammars/L1_trans_strands_4transpose.eg evidence_L1_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L1_trans_strands.stk -ar > garbage											#ancestral reconstruction
date >> times_L1_trans_strands																			#end time

date >> times_L1utr_trans_strands																				#start time
xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg evidence_L1utr_trans_strands.stk -t trained.eg	#train the grammar parameters
xrate -g trained.eg evidence_L1utr_trans_strands.stk -ar > garbage												#ancestral reconstruction
date >> times_L1utr_trans_strands																				#end time

done

