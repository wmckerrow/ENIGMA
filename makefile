all:
	make eix2stadn
	make GeneratePositions
	make L0_2_L1
	make L1_2_stadn
	make L1stk2L0stk
	make L1stk2L1UTRstk
	make L1utr_2_stadn
	make newL0stk2oldL0stk
	make oldL0stk2newL0stk
	make stk2gff_strands
	make strandrandomly
	make strandrandomly_L1
	make strandrandomly_L1utr
	make Generator
	make GeneratorL1
	make GeneratorL1_constgenes
	make Generator_L1UTR
	make maketranspose
	make accuracy_L0
	make accuracy_L1
	make exonpps_L0
	make exonpps_L1
	make exactpps_L0
	make exactpps_L1
	make simplify_ar

eix2stadn: converters/eix2stadn.cpp
	g++ converters/eix2stadn.cpp -o bin/eix2stadn

GeneratePositions: converters/GeneratePositions.cpp
	g++ converters/GeneratePositions.cpp -o bin/GeneratePositions

L0_2_L1: converters/L0_2_L1.cpp
	g++ converters/L0_2_L1.cpp -o bin/L0_2_L1

L1_2_stadn: converters/L1_2_stadn.cpp
	g++ converters/L1_2_stadn.cpp -o bin/L1_2_stadn

L1stk2L0stk: converters/L1stk2L0stk.cpp
	g++ converters/L1stk2L0stk.cpp -o bin/L1stk2L0stk

L1stk2L1UTRstk: converters/L1stk2L1UTRstk.cpp
	g++ converters/L1stk2L1UTRstk.cpp -o bin/L1stk2L1UTRstk

L1utr_2_stadn: converters/L1utr_2_stadn.cpp
	g++ converters/L1utr_2_stadn.cpp -o bin/L1utr_2_stadn

newL0stk2oldL0stk: converters/newL0stk2oldL0stk.cpp
	g++ converters/newL0stk2oldL0stk.cpp -o bin/newL0stk2oldL0stk

oldL0stk2newL0stk: converters/oldL0stk2newL0stk.cpp
	g++ converters/oldL0stk2newL0stk.cpp -o bin/oldL0stk2newL0stk

stk2gff_strands: converters/stk2gff_strands.cpp
	g++ converters/stk2gff_strands.cpp -o bin/stk2gff_strands

strandrandomly: converters/strandrandomly.cpp
	g++ converters/strandrandomly.cpp -o bin/strandrandomly

strandrandomly_L1: converters/strandrandomly_L1.cpp
	g++ converters/strandrandomly_L1.cpp -o bin/strandrandomly_L1

strandrandomly_L1utr: converters/strandrandomly_L1utr.cpp
	g++ converters/strandrandomly_L1utr.cpp -o bin/strandrandomly_L1utr

Generator: generators/Generator.cpp
	g++ generators/Generator.cpp -o bin/Generator

GeneratorL1: generators/GeneratorL1.cpp
	g++ generators/GeneratorL1.cpp -o bin/GeneratorL1

GeneratorL1_constgenes: generators/GeneratorL1_constgenes.cpp
	g++ generators/GeneratorL1_constgenes.cpp -o bin/GeneratorL1_constgenes

Generator_L1UTR: generators/Generator_L1UTR.cpp
	g++ generators/Generator_L1UTR.cpp -o bin/Generator_L1UTR

maketranspose: grammars/maketranspose.cpp
	g++ grammars/maketranspose.cpp -o bin/maketranspose

accuracy_L0: scoring_xrateout/accuracy_L0.cpp
	g++ scoring_xrateout/accuracy_L0.cpp -o bin/accuracy_L0

accuracy_L1: scoring_xrateout/accuracy_L1.cpp
	g++ scoring_xrateout/accuracy_L1.cpp -o bin/accuracy_L1

exonpps_L0: scoring_xrateout/exonpps_L0.cpp
	g++ scoring_xrateout/exonpps_L0.cpp -o bin/exonpps_L0

exonpps_L1: scoring_xrateout/exonpps_L1.cpp
	g++ scoring_xrateout/exonpps_L1.cpp -o bin/exonpps_L1

exactpps_L0: scoring_xrateout/exactpps_L0.cpp
	g++ scoring_xrateout/exactpps_L0.cpp -o bin/exactpps_L0

exactpps_L1: scoring_xrateout/exactpps_L1.cpp
	g++ scoring_xrateout/exactpps_L1.cpp -o bin/exactpps_L1

simplify_ar: scoring_xrateout/simplify_ar.cpp
	g++ scoring_xrateout/simplify_ar.cpp -o bin/simplify_ar

clean:
	mv bin/README temp1234567890
	rm bin/*
	mv temp1234567890 bin/README
