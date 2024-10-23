#-lpthread
g++ ./cpp/main.cpp -std=c++17 -O3 -o MAPSE

ulimit -s unlimited

# "ulimit -s unlimited" for generating queries on the dense dataset Prosper


#### k =0: Query Generation
#./MAPSE query/ctu network/ctu 0

#./MAPSE query/prosper network/prosper 0

#./MAPSE query/btc_2011 network/btc_2011 0

#./MAPSE query/dot network/dot 0

#./MAPSE query/terraforms network/terraforms 0

#./MAPSE query/bayc network/bayc 0

#./MAPSE query/meebits network/meebits 0

#./MAPSE query/artblock network/artblock 0




##### k = 1: Overall                   
#./MAPSE query/ctu network/ctu 1

#./MAPSE query/prosper network/prosper 1

#./MAPSE query/btc_2011 network/btc_2011 1

./MAPSE query/bayc network/bayc 1


##### k = 2: Varying Delta 


#./MAPSE query/ctu network/ctu 2

#./MAPSE query/prosper network/prosper 2

#./MAPSE query/btc_2011 network/btc_2011 2

#./MAPSE query/bayc network/bayc 2


##### k = 3: Investigating Num of Incremental Cases

#./MAPSE query/ctu network/ctu 3

#./MAPSE query/prosper network/prosper 3

#./MAPSE query/bayc network/bayc 3


##### k = 4: TransTime                   
#./MAPSE query/ctu network/ctu 4

#./MAPSE query/prosper network/prosper 4

#./MAPSE query/btc_2011 network/btc_2011 4

#./MAPSE query/bayc network/bayc 4


##### k = 5: Varying Size                   
#./MAPSE query/ctu network/ctu 5

#./MAPSE query/prosper network/prosper 5

#./MAPSE query/btc_2011 network/btc_2011 5
