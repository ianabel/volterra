
Volterra: Volterra.hpp main.cpp Volterra.cpp
	g++ -std=c++20 -g -O0 -I /usr/include/eigen3 -o Volterra main.cpp Volterra.cpp

gk_volterra: Volterra.hpp Volterra.cpp GK.hpp gk_volterra.cpp
	g++ -std=c++20 -g -O0  -I /usr/include/eigen3 -o gk_volterra gk_volterra.cpp Volterra.cpp

gamma_check: Volterra.hpp Volterra.cpp GK.hpp gamma_check.cpp
	g++ -std=c++20  -g -O0  -I /usr/include/eigen3 -o gamma_check gamma_check.cpp Volterra.cpp
