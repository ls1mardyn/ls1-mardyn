P ?= 14
N = $(shell perl -e 'print (2**$(P))')

$(info P = $(P))
$(info N = $(N))

.PHONY:default all

default: inputfiles

all: animake mktcts mkesfera mkcp
inputfiles: animake_input mkTcTS_input mkesfera_input mkCP_input

animake:
	make -C ../../tools/animake/
mktcts:
	make -C ../../tools/mktcts/
mkesfera:
	make -C ../../tools/mkesfera
mkcp:
	make -C ../../tools/mkCP


animake_input: animake
	../../tools/animake/animake M$(P)-animake -c 0.0015046 -f EOX -N $(N) -T 0.0011876 -u
mkTcTS_input: mktcts
	../../tools/mktcts/mkTcTS M$(P)-mkTcTS -c 0.06482 -d 0.6223 -N $(N) -R 2.5 -T 0.95 -u
mkCP_input: mkcp
	../../tools/mkCP/mkCP -P M$(P)-mkCP -a 0.001 -c 0.1 -d 3 -f Ar -h 125 -H 1 -LN $(shell perl -e 'print int($(N)/2)') -T 1 -uU 0 -V 0.25 -3 0.1
#../../tools/mkCP/mkCP -P M26-mkCP -a 0.001 -c 0.1 -d 3 -f Ar -h 125 -H 1 -LN 37400000 -T 1 -uU 0 -V 0.25 -3 0.1
mkesfera_input: mkesfera
	../../tools/mkesfera/mkesfera M$(P)-mkesfera -I 234.3 -i 0.6231 -O $(shell perl -e 'print ($(N)**(1/3) / 0.98087)') -o 0.06516 -T 0.95 -u
#../../tools/mkesfera/mkesfera M26-mkesfera -I 234.3 -i 0.6231 -O 414.3 -o 0.06516 -T 0.95 -u
