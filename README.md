analyzeskills
=============

Simple tool that decomposes response patterns into a skill set
and a positive disjunctive normal form-matrix that expresses 
which skill sets solve which item.


Prerequisites:
--------------
   libCPROPS 0.1.12,
   use archive found in ./libs.
   Build and install first.

Build:
------

```
 mkdir build
 cd build
 cmake ..
 make
```

Demo:
-----
```
 ./analyzeskills ../summerschool.txt
```

Output:
```
  >   Response Dimensions: 9 Subjects x 8 Items
  >   M=
  >   ........
  >   X.......
  >   ..XXX.XX
  >   ....XXX.
  >   .XXXX...
  >   .XXXXXX.
  >   .XX..X..
  >   XXXXXXX.
  >   XXXXXXXX
  >   Need 7 incomparable skill vectors.
  >   Relation Product Factorization A*B=M: 5 Factors
  >   A=
  >   ....X
  >   ..XX.
  >   .XXX.
  >   .X.X.
  >   XX.X.
  >   X.X..
  >   XX...
  >   .X...
  >   B=
  >   .....
  >   ....X
  >   .X...
  >   X....
  >   ...X.
  >   X.XX.
  >   ..X..
  >   X.XXX
  >   XXXXX
  >   Searching for pDNF L1 and skill matrix S1 with 4 elements...
  >   Largest antichain in 4-elementary powerset: 10
  >   L1( 0,.) = ..XX
  >   L1( 1,.) = X..X
  >           \/ X.X.
  >           \/ XX..
  >           \/ .XXX
  >   L1( 2,.) = .X.X
  >           \/ X..X
  >           \/ X.X.
  >           \/ XX..
  >   L1( 3,.) = .X.X
  >           \/ X..X
  >           \/ X.X.
  >   L1( 4,.) = .X.X
  >           \/ .XX.
  >           \/ X..X
  >           \/ X.X.
  >   L1( 5,.) = .XX.
  >           \/ X.X.
  >           \/ XX..
  >   L1( 6,.) = .X.X
  >           \/ .XX.
  >           \/ X.X.
  >   L1( 7,.) = .X.X
  >   S1=
  >   ....
  >   ..XX
  >   .X.X
  >   .XX.
  >   X..X
  >   X.X.
  >   XX..
  >   X.XX
  >   .XXX
  >   Searching for pDNF L2 and skill matrix S2 with 3 elements...
  >   Largest antichain in 3-elementary powerset: 3
  >   Too many incomparable responses (7) for 3 skills!
```

Response Pattern Format:
------------------------
  Each line gives the response of one subject,
  write 'X' for items that have been solved,
  otherwise, write '.'. No extra characters 
  allowed!



 ~~ Immo.
