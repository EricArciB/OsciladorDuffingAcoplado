# === Makefile ===
CXX       = g++
CXXFLAGS  = -std=c++17 -Wall
INCLUDES  = -Iinclude
SRC       = src/main.cpp src/rungekuta.cpp src/utilidad.cpp
TARGET    = bin/ocilador

# AUXILIAR
AUX_SRC   = scripts/map.cpp src/rungekuta.cpp
AUX_TARGET= bin/auxiliar

all: $(TARGET) $(AUX_TARGET)

$(TARGET): $(SRC)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

$(AUX_TARGET): $(AUX_SRC)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(AUX_TARGET) results/*.dat results/*.png maps/*.png

view:
	start maps\*.png
