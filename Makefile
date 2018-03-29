#Makefile for cmlsynth

lib.name = cmlsynth


cmlsynth~.class.sources = cmlsynth.c 

datafiles =		 README.txt \
			 LICENSE.txt

suppress-wunused = yes

include Makefile.pdlibbuilder
