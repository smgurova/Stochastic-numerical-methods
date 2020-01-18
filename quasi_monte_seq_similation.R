# Quasi - Monte Carlo method
# https://en.wikipedia.org/wiki/Quasi-Monte_Carlo_method

# insert some packages for ggplot

#generate a 10 point van der Corput sequence.
npoints = 10
qmc <- data.frame(j = 1:npoints, x = halton(npoints))
ggplot(qmc,aes(x=x,y=j)) + geom_point(size=4) + scale_y_continuous('Sample Number')
ggplot(qmc,aes(x=j,y=x)) + geom_line() + geom_point(size=4) + scale_x_continuous('Sample Number')

#Halton Sequence in 2-D
npoints = 50
xy = halton(npoints,dim=2)
qmc <- data.frame(j = 1:npoints)
qmc[,c("x","y")] = xy
ggplot(qmc,aes(x=x,y=y)) + geom_point(size=4)

#npoints = 50
xy = halton(npoints,dim=6)
qmc <- data.frame(j = 1:npoints)
qmc[,2:7] = xy
ggpairs(qmc,columns=2:7)Halton Sequence in 6-D

#2-D the Sobol sequence is similar to Halton.
npoints = 50
xy = sobol(npoints,dim=2)
qmc <- data.frame(j = 1:npoints)
qmc[,c("x","y")] = xy
ggplot(qmc,aes(x=x,y=y)) + geom_point(size=4)

#In 1-D the torus is the fractional part of :
npoints = 10
qmc <- data.frame(j = 1:npoints, x = torus(npoints), y = (1:npoints)*sqrt(2))
ggplot(qmc,aes(x=x,y=y,color=j)) + geom_point(size=4) + scale_y_continuous('j times root 2') + scale_colour_gradient(low="blue",high="red")

# 2D Torus seq
npoints = 50
xy = torus(npoints,dim=2)
qmc <- data.frame(j = 1:npoints)
qmc[,c("x","y")] = xy
ggplot(qmc,aes(x=x,y=y)) + geom_point(size=4)
