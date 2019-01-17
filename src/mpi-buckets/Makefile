all: tests
clean:
	rm src/*.gch
lib: src/*.h
	cd src && for header in *.h; do c++ -I/usr/include/mpich -std=c++11 -c -O3 $$header; done
tests: lib ./tests/*.cc
	cd tests &&                    \
	for f in *.cc;	               \
	do                             \
		d=`echo $$f | sed s/.cc/.out/`;	      \
		mpic++ -I../src -std=c++17 -O3 $$f -o ../build/$$d;    \
	done
