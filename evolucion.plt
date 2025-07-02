set term gif animate
set output 'animacion_2.gif'

set xra[-10:1000]
set yra[0:2]

#set view map

set arrow 1 from 400,0 to 400,1 lc rgb '#00000' nohead
set arrow 3 from 400,1 to 600,1 lc rgb '#00000' nohead
set arrow 2 from 600,1 to 600,0 lc rgb '#00000' nohead
set arrow 4 from 0,0 to 0,2 lc rgb '#00000' nohead
set arrow 5 from 1000,0 to 1000,2 lc rgb '#00000' nohead

do for[a=1:10000]{
   plot 'funcion_onda.dat' index a u 1:2 w l t '',0.0 t '' lc rgb '#00000'
}

#set terminal png size 1024,720 enhanced font "Helvetica,15"
#set output 'escalon.png'


#plot 'funcion_onda.dat' i 0 u 1:2 w l t 'Partícula'
