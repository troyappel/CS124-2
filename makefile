strassen: main.cc
	g++ -o strassen main.cc -O3 -pthread

run: clean strassen

clean:
	rm -f randmst
