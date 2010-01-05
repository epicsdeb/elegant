# clean all the output files 
clean:
	rm -f Pelegant.* elegant.*  *.w1 *.w2 *.out *.bun
# compare the parallel version and serial version of the elegant
compare:
	sddsdiff manyParticles_p.w1 manyParticles_s.w1
	sddsdiff manyParticles_p.w2 manyParticles_s.w2
	sddsdiff manyParticles_p.out manyParticles_s.out
