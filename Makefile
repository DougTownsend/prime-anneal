
test_eq: src/eq.cpp src/test_eq.cpp include/eq.hpp src/thal.cpp include/thal.h
	g++ -Wall src/eq.cpp src/test_eq.cpp src/thal.cpp -Iinclude -L/share/tuck/dktownse/usr/lib -I/share/tuck/dktownse/usr/include -O3 -lm -lpthread -lmpfr -lgmp -mavx -o test_eq

clean:
	rm test_eq