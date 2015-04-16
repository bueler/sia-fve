
mahaffy:
	ln -s ../mahaffy

# for generating dome verification convergence figure in ../paper/
dome.verif: mahaffy domeconv.sh
	./domeconv.sh dome.verif 6 6

true.verif: mahaffy domeconv.sh
	./domeconv.sh true.verif 4 6 "-mah_true -snes_fd_color -da_grid_x 20 -da_grid_y 20"

# for generating bedstep verification table in ../paper/
bedstep.verif: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstep.verif 4 6 "-snes_max_it 200 -cs_end 12"

bedstepnoup.verif: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstepnoup.verif 4 6 "-snes_max_it 200 -cs_end 12 -mah_lambda 0.0"

bedstepfullup.verif: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstepfullup.verif 4 6 "-snes_max_it 200 -cs_end 12 -mah_lambda 1.0"

