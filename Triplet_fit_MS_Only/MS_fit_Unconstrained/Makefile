#
# makefile
#
# 
#

# files
srcs =triplet_MSfit_Unconst.cc
exe = $(srcs:.cc=)

# libs

cxxflags += -g -Os -Wall -Wno-deprecated -fPIC -pipe
cxxflags += $(shell root-config --cflags)

ldflags  += $(shell root-config --ldflags --libs)

# targets

all	: start $(exe) end

start	:
	@echo ''
	@echo ' * make ... '
	@echo ''
	@rm -f ./#* ./*~ ./*.*~

%	: %.cc	
	$(CXX) $(cxxflags) -o $@ $< $(ldflags)

end	:
	@echo ' * Good'
	@echo ''

clean	:
	@echo ''
	@echo ' * Next Step'
	@echo ''
	$(RM) $(exe)
