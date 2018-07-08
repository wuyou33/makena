.DEFAULT_GOAL := all

# Change the following two locations to yours.
GOOGLE_TEST_INC_DIR   = -I./GoogleTest/include
GOOGLE_TEST_LIB_DIR   = -L./GoogleTest/lib

INSTALL_DIR_LIB = /usr/local/lib
INSTALL_DIR_BIN = /usr/local/bin

.PHONY:	unit_tests
.PHONY: interactive_tests
.PHONY: all
.PHONY: install
.PHONY: clean

CC           = c++
LD           = c++
RM           = rm -f
RMR          = rm -fr
CP           = cp
PYTHON       = python

PWD          =  $(shell pwd)
SRC_DIR_LIB  = ./src_lib
SRC_DIR_UT   = ./src_unit_tests
SRC_DIR_IT   = ./src_interactive_tests
SRC_DIR_BIN  = ./src_bin
BIN_DIR      = ./bin
LIB_DIR      = ./libs

OBJ_DIR_REL  = ./objs_release
OBJ_DIR_DBG  = ./objs_debug
OBJ_DIR_UT   = ./objs_unit_tests
BIN_DIR_IT   = ./bins_interactive_tests

CPPFLAGS     = -Wall -std=c++1y
CPPFLAGS_REL = -fPIC -O3
CPPFLAGS_DBG = -g -O0 -DUNIT_TESTS
CPPFLAGS_BIN = -L./libs

CFLAGS_DBG   = $(GOOGLE_TEST_INC_DIR) -I. -I./include -I/usr/local/include/wailea -I/usr/local/include/makena
CFLAGS_REL   = -I. -I./include -I/usr/local/include/wailea -I/usr/local/include/makena

LDFLAGS      = -shared -lwailea
TARGET_LIB   = $(LIB_DIR)/libmakena.so

UNIT_TESTER  = unit_tester

INTERACTIVE_TEST_LIBS = -lwailea -framework CoCoa -framework OpenGL -framework IOKit -framework CoreFoundation -framework CoreVideo -lAntTweakBar -lglfw3 -lGLEW -lstdc++

DIR_GUARD    = @mkdir -p $(@D)

OBJS_REL     = $(patsubst %,$(OBJ_DIR_REL)/%, \
                         $(subst .cpp,.o, \
                           $(notdir \
                             $(wildcard \
                               $(SRC_DIR_LIB)/*.cpp ))))

OBJS_DBG     = $(patsubst %,$(OBJ_DIR_DBG)/%, \
                         $(subst .cpp,.o, \
                           $(notdir \
                             $(wildcard \
                               $(SRC_DIR_LIB)/*.cpp ))))

OBJS_UT      = $(patsubst %,$(OBJ_DIR_UT)/%, \
                         $(subst .cpp,.o, \
                           $(notdir \
                             $(wildcard \
                               $(SRC_DIR_UT)/*.cpp ))))

BINS_IT      = $(patsubst %,$(BIN_DIR_IT)/%, \
                         $(basename \
                           $(notdir \
                             $(wildcard \
                               $(SRC_DIR_IT)/*.cpp ))))


$(OBJ_DIR_REL)/%.o: $(SRC_DIR_LIB)/%.cpp
	$(DIR_GUARD)
	$(CC) $(CFLAGS_REL) $(CPPFLAGS) $(CPPFLAGS_REL) -c $< -o $@

$(OBJ_DIR_DBG)/%.o: $(SRC_DIR_LIB)/%.cpp
	$(DIR_GUARD)
	$(CC) $(CFLAGS_DBG) $(CPPFLAGS) $(CPPFLAGS_DBG) -c $< -o $@

$(OBJ_DIR_UT)/%.o: $(SRC_DIR_UT)/%.cpp
	$(DIR_GUARD)
	$(CC) $(CFLAGS_DBG) $(CPPFLAGS) $(CPPFLAGS_DBG) -c $< -o $@

$(TARGET_LIB):	$(OBJS_REL)
	$(DIR_GUARD)
	$(LD) $(LDFLAGS) -o $@ $^

$(UNIT_TESTER): $(OBJS_DBG) $(OBJS_UT)
	$(LD) $(GOOGLE_TEST_LIB_DIR) $^ -lwailea -lgtest -lgtest_main -lstdc++ -o $@

$(BIN_DIR_IT)/%:	$(SRC_DIR_IT)/%.cpp $(OBJS_DBG)
	$(DIR_GUARD)
	$(CC) $(CFLAGS_DBG) $(CPPFLAGS) $(CPPFLAGS_DBG) $^ $(INTERACTIVE_TEST_LIBS) -o $@

unit_tests:	$(UNIT_TESTER)
	./$^

interactive_tests:	$(BINS_IT)

all:	$(TARGET_LIB)

clean:
	-$(RM) $(TARGET_LIB) $(UNIT_TESTER)
	-$(RMR) $(OBJ_DIR_REL) $(OBJ_DIR_DBG) $(OBJ_DIR_UT) $(LIB_DIR) $(BIN_DIR_IT)
	-$(RM) include/gjk_origin_finder_auto_generated.hpp

install:
	sudo $(CP) $(TARGET_LIB) $(INSTALL_DIR_LIB)
