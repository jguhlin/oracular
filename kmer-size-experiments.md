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