#-lpthread
g++ ./cpp/main.cpp -std=c++17 -O3 -o BFA

# "ulimit -s unlimited" for generating queries on the dense dataset Prosper


##### k = 1: Overall                   
#./BFA query/ctu network/ctu 1

#./BFA query/prosper network/prosper 1

#./BFA query/btc_2011 network/btc_2011 1




##### k = 2: Varying Delta 


#./BFA query/ctu network/ctu 2

#./BFA query/prosper network/prosper 2

#./BFA query/btc_2011 network/btc_2011 2




##### k = 3: Investigating Num of Incremental Cases

#./BFA query/ctu network/ctu 3

#./BFA query/prosper network/prosper 3




##### k = 4: TransTime                   
#./BFA query/ctu network/ctu 4

#./BFA query/prosper network/prosper 4

#./BFA query/btc_2011 network/btc_2011 4




##### k = 5: Varying Size                   
#./BFA query/ctu network/ctu 5

#./BFA query/prosper network/prosper 5

#./BFA query/btc_2011 network/btc_2011 5
