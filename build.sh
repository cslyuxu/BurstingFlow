#-lpthread
g++ ./cpp/main.cpp -std=c++17 -O3 -o BFQ

ulimit -s unlimited

# "ulimit -s unlimited" for generating queries on the dense dataset Prosper


#### k =0: Query Generation
#./BFQ query/ctu network/ctu 0

#./BFQ query/prosper network/prosper 0

#./BFQ query/btc_2011 network/btc_2011 0

#./BFQ query/dot network/dot 0

#./BFQ query/terraforms network/terraforms 0

#./BFQ query/bayc network/bayc 0

#./BFQ query/meebits network/meebits 0

#./BFQ query/artblock network/artblock 0




##### k = 1: Overall                   
#./BFQ query/ctu network/ctu 1

#./BFQ query/prosper network/prosper 1

#./BFQ query/btc_2011 network/btc_2011 1

./BFQ query/bayc network/bayc 1


##### k = 2: Varying Delta 


#./BFQ query/ctu network/ctu 2

#./BFQ query/prosper network/prosper 2

#./BFQ query/btc_2011 network/btc_2011 2

#./BFQ query/bayc network/bayc 2


##### k = 3: Investigating Num of Incremental Cases

#./BFQ query/ctu network/ctu 3

#./BFQ query/prosper network/prosper 3

#./BFQ query/bayc network/bayc 3


##### k = 4: TransTime                   
#./BFQ query/ctu network/ctu 4

#./BFQ query/prosper network/prosper 4

#./BFQ query/btc_2011 network/btc_2011 4

#./BFQ query/bayc network/bayc 4


##### k = 5: Varying Size                   
#./BFQ query/ctu network/ctu 5

#./BFQ query/prosper network/prosper 5

#./BFQ query/btc_2011 network/btc_2011 5
