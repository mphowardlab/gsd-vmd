display resetview
color Display Background white
axes location Off

mol delrep 0 0

mol color Name
mol representation VDW 1.000000 12.000000
mol material AOChalky
mol selection all
mol addrep 0
mol representation Bonds 0.300000 12.000000
mol addrep 0

pbc box -width 1 -color black -center origin

scale by 0.5
translate by 0.1 0 0
