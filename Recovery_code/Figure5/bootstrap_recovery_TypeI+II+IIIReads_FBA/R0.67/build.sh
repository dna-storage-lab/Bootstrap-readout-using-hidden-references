
# ----------------------------------------------------------
# Compile Programs
# ----------------------------------------------------------
cd ./src/
gcc -o ../bin/Seq_chunker Seq_chunker.c
gcc -o ../bin/SubsampleFastqRandom SubsampleFastqRandom.c
gcc -o ../bin/CalBitError CalBitError.c

# ----------------------------------------------------------
# Compile step1 Programs
# ----------------------------------------------------------
gcc -o ../bin/SlidingCorrelation SlidingCorrelation.c -lpthread
gcc -o ../bin/getthre getthre.c
gcc -o ../bin/post_dec_hamming_dis post_dec_hamming_dis.c

# ----------------------------------------------------------
# Compile step2 Programs
# ----------------------------------------------------------
gcc -o ../bin/majorityvoting majorityvoting.c
c++ -c ../lib/edlib/src/edlib.cpp -o ../bin/edlib.o -I ../lib/edlib/include
cc -c lowthres_pthread_edlib.c -o ../bin/lowthres_pthread_edlib.o -I ../lib/edlib/include -lpthread
c++ ../bin/lowthres_pthread_edlib.o ../bin/edlib.o -o ../bin/lowthres_pthread_edlib -lpthread

# ----------------------------------------------------------
# Compile step3 Programs
# ----------------------------------------------------------
g++ -o ../bin/decode_feedback_align  decode_feedback_align.cpp ../lib/edlib/src/edlib.cpp -I ../lib/edlib/include