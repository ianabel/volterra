
Volterra: Volterra.hpp main.cpp Volterra.cpp
	g++ -I /usr/include/eigen3 -o Volterra main.cpp Volterra.cpp

gk_volterra: Volterra.hpp Volterra.cpp GK.hpp gk_volterra.cpp
	g++ -I /usr/include/eigen3 -o gk_volterra gk_volterra.cpp Volterra.cpp
