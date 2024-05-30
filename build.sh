#-lpthread
g++ ./cpp/main.cpp -std=c++17 -O3 -o BFQ

# "ulimit -s unlimited" for generating queries on the dense dataset Prosper


##### k = 1: Overall                   
#./BFQ query/ctu network/ctu 1

#./BFQ query/prosper network/prosper 1

#./BFQ query/btc_2011 network/btc_2011 1




##### k = 2: Varying Delta 


#./BFQ query/ctu network/ctu 2

#./BFQ query/prosper network/prosper 2

#./BFQ query/btc_2011 network/btc_2011 2




##### k = 3: Investigating Num of Incremental Cases

#./BFQ query/ctu network/ctu 3

#./BFQ query/prosper network/prosper 3




##### k = 4: TransTime                   
#./BFQ query/ctu network/ctu 4

#./BFQ query/prosper network/prosper 4

#./BFQ query/btc_2011 network/btc_2011 4




##### k = 5: Varying Size                   
#./BFA query/ctu network/ctu 5

#./BFA query/prosper network/prosper 5

#./BFA query/btc_2011 network/btc_2011 5
