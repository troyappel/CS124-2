strassen: strassen.cpp
	g++ -std=c++17 -O2 -Wall -Wextra strassen.cpp -o strassen -lm -lpthread

run: clean strassen exec

clean:
	rm -f strassen

# define make_input = 
# 	#!/bin/bash
# 	> input.txt ; 
# 	awk -v n=2000000 'BEGIN{srand(); for(i=1;i<=n;i++) print(int(rand()*10))}' > input.txt
# endef



exec:
	# $(make_input) ;
	# num1=2; \
	# while [ "$$num1" -le 1024 ] ; do \
	# 	echo $$num1 ; \
	# 	./strassen $$num1 1024 input.txt ; \
	# 	num1=$$((num1 * 2)) ; \
	# done ; \
	# true
	# ./strassen 64 1024 input.txt
.ONESHELL: