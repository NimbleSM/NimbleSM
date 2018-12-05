all: tests
lib: src/*.hpp
	echo "Made lib"
tests: lib ./tests/*.cc
	cd tests &&                    \
	for f in *.cc;	               \
	do                             \
		d=`echo $$f | sed s/.cc/.out/`;	      \
		mpic++ -I../src -std=c++11 -O3 $$f -o ../build/$$d;    \
	done
