include_directories(${BICEPS_SOURCE_DIR}/include ${BICEPS_BINARY_DIR}/include)


SET(PEPSPLICE_SRCS
	BestMatches.cpp
	Chromosome.cpp
	Chromosomes.cpp
	Distribution.cpp
	DnaAA.cpp
	HotSpectra.cpp
	HotSpots.cpp
	Hypergeometric.cpp
	Match.cpp
	Peak.cpp
	Protein.cpp
	ProteinParser.cpp
	Results.cpp
	ResultsPattern.cpp
	Scoring.cpp
	Services.cpp
	SlidingWindow.cpp
	Spectra.cpp
	Spectrum.cpp
	SpectrumParser.cpp
	Timer.cpp
	Tuple.cpp
	Tuples.cpp
	UnreachableBins.cpp
    pepsplice.cpp
)
add_library(pepsplice STATIC ${PEPSPLICE_SRCS})
