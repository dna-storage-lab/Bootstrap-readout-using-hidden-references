
# ----------------------------------------------------------
# Compile Programs
# ----------------------------------------------------------
cd ./src/

gcc -o ../bin/Seq_chunker Seq_chunker.c
gcc -o ../bin/SlidingCorrelation SlidingCorrelation.c -lpthread
gcc -o ../bin/BitwiseConsensusRecovery BitwiseConsensusRecovery.c
gcc -o ../bin/PostDecHammingDistance PostDecHammingDistance.c
gcc -o ../bin/SubsampleFastqRandom SubsampleFastqRandom.c
gcc -o ../bin/CalBitError CalBitError.c
