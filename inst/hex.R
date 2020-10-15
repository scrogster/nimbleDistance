#code to make a package hexsticker
library(hexSticker)

plotfunc<-function(){
	plot(y=dnorm(seq(0, 5, by=0.01), 0, 1), x=seq(0, 5, by=0.01), cex.lab=1.5, las=1,
			 type="l", xlab="x", ylab="p(x)", col="red", lwd=2, mai=c(2,0,0,0),  xaxt='n', yaxt="n", ann=FALSE)
	lines(x=rnorm(500, 1, 0.22), y=seq(0, dnorm(0, 0,1), length.out = 500), col="grey", lwd=0.4)
	abline(v=1, col="red", lty=2)
	rug(abs(rnorm(50, 0, 1)), col="wheat3", ticksize=0.005)
	mtext(expression(sigma), 1, at=1, cex=3, line=-0.3, col="wheat3")
	box(col="wheat3", lwd=0.4)
}

s <- sticker(~plotfunc(),
						 package="nimbleDistance", p_size=18, p_color ="wheat3",
						 s_x=.8, s_y=.6, s_width=1.6, s_height=1.6,
						 h_fill="#331188", h_color = "wheat3", h_size=0.9,
						 filename="inst/nimbleDistance.png")



