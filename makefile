strassen: main.cc
	g++ -o strassen main.cc -O3 -pthread

run: clean strassen exec

clean:
	rm -f randmst

exec:
	./strassen