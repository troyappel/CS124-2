strassen: main.cc
	g++ -o strassen main.cc -O1 -pthread -g

run: clean strassen exec

clean:
	rm -f strassen

define make_input = 
	#!/bin/bash
	> input.txt ; 
	awk -v n=1000 'BEGIN{srand(); for(i=1;i<=n;i++) print(int(rand()*10))}' > input.txt
endef

exec:
	$(make_input) ; 
	./strassen 0 15 input.txt

.ONESHELL: