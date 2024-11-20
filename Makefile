
test_eq: src/eq.cpp src/test_eq.cpp include/eq.hpp src/thal.cpp include/thal.h
	g++ -Wall src/eq.cpp src/test_eq.cpp src/thal.cpp -Iinclude -O3 -lm -lpthread -lmpfr -lgmp -march=native -o test_eq

clean:
	rm eq_ratio calc_eq test_eq