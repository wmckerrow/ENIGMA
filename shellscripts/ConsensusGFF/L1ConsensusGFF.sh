rm exon_sens_spec

dotplotter $ENIGMADIR/RealData/newgenomes_0423.maf $ENIGMADIR/RealData/info_files/ArabFull.info Arabidopsis $ENIGMADIR/RealData/info_files/MediFull.info Medicago $ENIGMADIR/RealData/info_files/SoyFull.info Soybean Arabidopsis Medicago > AM_raw.dot
getdotsubset AM_raw.dot 1 5000000 1 5000000 25 > AM_25.dot

dotplotter $ENIGMADIR/RealData/newgenomes_0423.maf $ENIGMADIR/RealData/info_files/ArabFull.info Arabidopsis $ENIGMADIR/RealData/info_files/MediFull.info Medicago $ENIGMADIR/RealData/info_files/SoyFull.info Soybean Soybean Medicago > SM_raw.dot
getdotsubset SM_raw.dot 1 5000000 1 5000000 25 > SM_25.dot

prepare4stk 1 A AM_25.dot 1 0 $ENIGMADIR/RealData/arabSynteny.gff M 4 $ENIGMADIR/RealData/medicago_evidence_synteny_est2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_protein2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_snap_masked.gff $ENIGMADIR/RealData/medicago_evidence_synteny_augustus_masked.gff 1 1500000 > arabSections
mv goodmafs0 AM.dot
prepare4stk 1 S SM_25.dot 1 0 $ENIGMADIR/RealData/soySynteny.gff M 4 $ENIGMADIR/RealData/medicago_evidence_synteny_est2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_protein2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_snap_masked.gff $ENIGMADIR/RealData/medicago_evidence_synteny_augustus_masked.gff 1 1500000 > soySections
mv goodmafs0 SM.dot
choosesections 2 arabSections soySections .1 > Sections

#medi current annotation for comparison by eval
getgffsections $ENIGMADIR/RealData/mediSynteny.gff all Sections > mediSections.gff
gff3_to_gtf.pl mediSections.gff > mediSections.gtf

#3way
sectionAlignmentsStk_L1 2 A AM.dot 1 0 $ENIGMADIR/RealData/arabSynteny.gff S SM.dot 1 0 $ENIGMADIR/RealData/soySynteny.gff M 4 $ENIGMADIR/RealData/medicago_evidence_synteny_est2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_protein2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_snap_masked.gff $ENIGMADIR/RealData/medicago_evidence_synteny_augustus_masked.gff L1.pos Sections "(M1:.1,M2:.1,M3:.1,M4:.1,((S1:.1)S:.5,((A1:.1)A:.9)AMS:.4)MS:.5)M" > step1.stk
L1stk2L1UTRstk step1.stk > step2.stk
evidencereadingframe step2.stk L1.pos M1 1 M2 2 M3 2 M4 2 > L1.stk
xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg L1.stk -x expanded.eg -gc "#"
maketranspose expanded.eg species reverse AMS MS > transpose.eg
xrate -g transpose.eg L1.stk -t trained.eg -gc "#"
xrate -g trained.eg L1.stk -gc "#" -ar > arout
simplify_ar arout > xrateout
stk2gff_L1utr 0 xrateout consensus3way L1.pos
gff3_to_gtf.pl consensus3way_M.gff > consensus3way_M.gtf
evaluate_gtf.pl mediSections.gtf consensus3way_M.gtf > evalout_3way
#echo "3way" >> exon_sens_spec
parseeval evalout_3way exon_sens_spec 21 2 196 3 509 3 518 3 526 4 790 3 794 3 > garbage

#2way
all2way SM_25.dot > SM_25_2way.dot
sectionAlignmentsStk_L1 1 S SM_2way.dot 1 0 $ENIGMADIR/RealData/soySynteny.gff M 4 $ENIGMADIR/RealData/medicago_evidence_synteny_est2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_protein2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_snap_masked.gff $ENIGMADIR/RealData/medicago_evidence_synteny_augustus_masked.gff L1.pos Sections "(M1:.1,M2:.1,M3:.1,M4:.1,((S1:.1)S:.5,((A1:.1)A:.9)AMS:.4)MS:.5)M" > step1.stk
L1stk2L1UTRstk step1.stk > step2.stk
evidencereadingframe step2.stk L1.pos M1 1 M2 2 M3 2 M4 2 > L1.stk
xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg L1.stk -x expanded.eg -gc "#"
maketranspose expanded.eg species reverse AMS MS > transpose.eg
xrate -g transpose.eg L1.stk -t trained.eg -gc "#"
xrate -g trained.eg L1.stk -gc "#" -ar > arout
simplify_ar arout > xrateout
stk2gff_L1utr 0 xrateout consensus2way L1.pos
gff3_to_gtf.pl consensus2way_M.gff > consensus2way_M.gtf
evaluate_gtf.pl mediSections.gtf consensus2way_M.gtf > evalout_2way
#echo "3way" >> exon_sens_spec
parseeval evalout_2way exon_sens_spec 21 2 196 3 509 3 518 3 526 4 790 3 794 3 > garbage

#1way
touch empty.dot
sectionAlignmentsStk_L1 1 S empty.dot 1 0 $ENIGMADIR/RealData/soySynteny.gff M 4 $ENIGMADIR/RealData/medicago_evidence_synteny_est2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_protein2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_snap_masked.gff $ENIGMADIR/RealData/medicago_evidence_synteny_augustus_masked.gff L1.pos Sections "(M1:.1,M2:.1,M3:.1,M4:.1,((S1:.1)S:.5,((A1:.1)A:.9)AMS:.4)MS:.5)M" > step1.stk
L1stk2L1UTRstk step1.stk > step2.stk
evidencereadingframe step2.stk L1.pos M1 1 M2 2 M3 2 M4 2 > L1.stk
xrate -g $ENIGMADIR/grammars/L1utr_trans_strands_4transpose.eg L1.stk -x expanded.eg -gc "#"
maketranspose expanded.eg species reverse AMS MS > transpose.eg
xrate -g transpose.eg L1.stk -t trained.eg -gc "#"
xrate -g trained.eg L1.stk -gc "#" -ar > arout
simplify_ar arout > xrateout
stk2gff_L1utr 0 xrateout consensus1way L1.pos
gff3_to_gtf.pl consensus1way_M.gff > consensus1way_M.gtf
evaluate_gtf.pl mediSections.gtf consensus2way_M.gtf > evalout_1way
#echo "3way" >> exon_sens_spec
parseeval evalout_1way exon_sens_spec 21 2 196 3 509 3 518 3 526 4 790 3 794 3 > garbage

gff3_to_gtf.pl consensus1way_M1.gff > consensus1way_M1.gtf
evaluate_gtf.pl mediSections.gtf consensus1way_M1.gtf > evalout_est
#echo "est stkline" >> exon_sens_spec
parseeval evalout_est exon_sens_spec 21 2 196 3 509 3 518 3 526 4 790 3 794 3 > garbage

gff3_to_gtf.pl consensus1way_M2.gff > consensus1way_M2.gtf
evaluate_gtf.pl mediSections.gtf consensus1way_M2.gtf > evalout_protein
#echo "protein stkline" >> exon_sens_spec
parseeval evalout_protein exon_sens_spec 21 2 196 3 509 3 518 3 526 4 790 3 794 3 > garbage

gff3_to_gtf.pl consensus1way_M3.gff > consensus1way_M3.gtf
evaluate_gtf.pl mediSections.gtf consensus1way_M3.gtf > evalout_snap
#echo "snap stkline" >> exon_sens_spec
parseeval evalout_snap exon_sens_spec 21 2 196 3 509 3 518 3 526 4 790 3 794 3 > garbage

gff3_to_gtf.pl consensus1way_M4.gff > consensus1way_M4.gtf
evaluate_gtf.pl mediSections.gtf consensus1way_M4.gtf > evalout_augustus
#echo "augustus stkline" >> exon_sens_spec
parseeval evalout_augustus exon_sens_spec 21 2 196 3 509 3 518 3 526 4 790 3 794 3 > garbage

scoreSectionDistances Sections 8 $ENIGMADIR/RealData/mediSynteny.gff consensus3way_M.gff consensus2way_M.gff consensus1way_M.gff $ENIGMADIR/RealData/medicago_evidence_synteny_est2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_protein2genome.gff $ENIGMADIR/RealData/medicago_evidence_synteny_snap_masked.gff $ENIGMADIR/RealData/medicago_evidence_synteny_augustus_masked.gff > scores
