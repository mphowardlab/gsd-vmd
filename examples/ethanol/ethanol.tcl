display resetview
color Display Background white
axes location Off

mol new ethanol.gsd type gsd first 0 last -1 step 1 waitfor 1
mol delrep 0 0

mol color Name
mol representation VDW 1.000000 12.000000
mol material AOChalky
mol selection all
mol addrep 0
mol representation Bonds 0.300000 12.000000
mol addrep 0

pbc box -width 1 -color black -center origin
