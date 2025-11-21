run: bin/Programa
	./bin/Programa

bin/Programa: src/Programa.cpp 
	@mkdir -p bin
	g++ -Iinclude -o bin/Programa src/Programa.cpp
