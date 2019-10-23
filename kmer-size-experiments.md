** Speed Experiments
at k=13, t1ha hash is the fastest, followed by wyhash

at k=11, t1ha is the fastest, followed by wyhash and fnv (tied, effectively)

Can use bitvec directly and get a u64 out of it (probably the best option) -- limits kmer sizes up to 21, though that's excessive anyway for our use..

** For NT
Kmer size of k=13, min_count = 10 RC NOT FIXED
35442141/118077428 == 30% of kmers used...
Min count of 50 =~ 10% of kmers used...
Avg. occurences at k=13 1715.6592

At k=11 , min_count = 10 RC FIXED
13086502/16898657 == 0.77441077122 == 77% of kmers have meaning
Avg. occurences at k=11 11994.443


** Kmer Size Experiments
Conclusion: k=13 for nt, k=11 for individual genome / reads

Kmer size experiments...

capitalize before good seq
77.03user 4.61system 1:21.76elapsed 99%CPU (0avgtext+0avgdata 19030228maxresident)k

capitalize after good seq
76.83user 4.61system 1:21.55elapsed 99%CPU (0avgtext+0avgdata 19030136maxresident)k


without RC check...
76.88user 4.76system 1:21.75elapsed 99%CPU (0avgtext+0avgdata 19030152maxresident)k

with RC check...
76.83user 4.61system 1:21.55elapsed 99%CPU (0avgtext+0avgdata 19030136maxresident)k


WITH rc and proper capitalization

k=13
    342981308
     52270360
     6.561679

k=12
    343080666
     16353228
    20.979385

k=11
    343180024
      4235398
    81.026634

k=21
    342186444
    319348610
    1.0715138

So the stats below (not timing, but kmer stats) are marred by an issue...

WITH reverse complement
   k=13 342981308
         66930112
         expected average: 5.12

   k=10 343279382 total
          5237227 unique
          expected average: 65.546

   k=11 343180024 total 
         11746932 unique
          expected average: 29.214

   k=9  343378740
          2235941
        expected average: 153.57


   These do not include reverse complement...
 * k=11 171590012 total tokens (This is what I suggest)
         10769668 unique kmers
         expected average: 15.93
   
   k=10 171639691
          4780515
         expected average: 35.9

   k=9  171689370 total
          2079595 unique
          expected avg: 82.559 (probably going towards useless?)

   k=25 171093222 total
        162723215 unique
        expected average: 1.05 (probably useless)

   k=17 171291938 total
        150872681 unique
        expected average: 1.135 (probably useless)

Speed tests:
Just straight add
k=13
153.88user 11.19system 2:45.25elapsed 99%CPU (0avgtext+0avgdata 12136604maxresident)k

// With RwLock but single-threaded...
166.41user 12.44system 2:59.09elapsed 99%CPU (0avgtext+0avgdata 18407048maxresident)k
0inputs+0outputs (0major+10285501minor)pagefaults 0swaps

// With RwLock but multi-threaded...

sort and add based on previous id !!! ... so no
9683.21user 686.98system 5:49.42elapsed 2967%CPU (0avgtext+0avgdata 12135664maxresident)k
0inputs+0outputs (0major+126610995minor)pagefaults 0swaps
